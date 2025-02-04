library(dplyr)
library(tidyr)
library(limma)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(stringr)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)

PATH_results = "./output/explants/"

df      <- read.csv("./data/explants/combined_genes.csv", row.names = 1, check.names = FALSE)
colData <- read.csv("./data/explants/colData.csv", row.names = 1)

enrichments <- read.csv("./output/explants/tissue-enr/enrichments.csv", row.names = 1)

#-------------------------------------------------------------------------------

head(enrichments)

# Extract neuronal enriched genes (see 06-explant-enrichment.R)
mat <- df[rownames(df) %in% enrichments$Gene, ]

colData <- colData[colData$sampleID %in% colnames(mat), ]
index   <- match(colnames(mat), colData$sampleID)
colData <- colData[index, ]

head(colnames(mat))
head(colData$sampleID)

TC <- colData$sampleID[colData$Turbo == "T"]

mat <- mat[, colnames(mat) %in% TC]
colData <- colData[colData$sampleID %in% TC, ]

design <- model.matrix(~ 0 + Condition, data = colData) #select design
head(design)

fit <- lmFit(mat, design)

colnames(fit) #check possible comparisons

#----------------

# Oxaliplatin
contrast_matrix <- makeContrasts(ConditionOx -  ConditionV, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust="BH", number=Inf)

head(results)

sig_de <- results[results$adj.P.Val < 0.05, , drop = FALSE] #extract significant proteins if there
head(sig_de)

sig_de <- na.omit(sig_de)
dim(sig_de)

# volcano plotting
res       <- as.data.frame(results) 
mutateddf <- mutate(res, Sig=ifelse(res$adj.P.Val<0.05 & abs(res$logFC)>0.5, "FDR < 0.05 & LFC > 0.5", ifelse("NS")))
input     <- cbind(gene=rownames(res), mutateddf) 

volc = ggplot(input, aes(logFC, -log10(P.Value))) + geom_point(aes(col=Sig)) +
  #scale_color_manual(values = c("#0D0887FF","#9512A1FF", "grey")) + 
  scale_color_manual(values = c("#9357f1", "grey")) + 
  ggrepel::geom_text_repel(data=subset(input, input$gene %in% rownames(sig_de)),
                           aes(label=gene), size=4, segment.alpha= 0.2, force =2, max.overlaps=16) 
volc <- volc + theme_bw() + theme(aspect.ratio=1)
volc <- volc + theme(legend.position="bottom", axis.text.y = element_text(size= 12, face="bold"), 
                     axis.title.y = element_text(size=14), axis.title.x = element_text(size= 14), 
                     axis.text.x = element_text(size= 12), legend.title=element_text(size=14), 
                     legend.text=element_text(size=14), plot.title=element_text(size=12, hjust = 0.5)) + 
  ggtitle("Oxaliplatin vs Vehicle")

print(volc)

pdf(file = paste0(PATH_results, "volcano.pdf"), height = 4, width = 4)
print(volc)
dev.off()

pdf(file = paste0(PATH_results, "volcano_big.pdf"), height = 6, width = 6)
print(volc)
dev.off()

table(sig_de$adj.P.Val < 0.05)

write.csv(results, paste(PATH_results, "DEP-analysis-limma_drugeffect.csv"))
write.csv(sig_de,  paste(PATH_results, "DEP-analysis-limma_drugeffect_sig.csv"))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Biological pathway differences
# NOTE: empty for UP

test_df    <- rownames(sig_de)[sig_de$logFC < -0.5] #adjust for up, down regulated genes
background <- rownames(results) # only test genes that could be tested with limma, not all enriched

ego <- enrichGO(gene          = test_df,
                universe      = background,
                OrgDb         = org.Mm.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

# remove redundancy in the GO terms
ego2 <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min, measure = 'Wang')

pdf(file = paste(PATH_results, "BP_UP-network.pdf", sep=""), width = 8, height = 8)
goplot(ego2)
dev.off()

pdf(file = paste(PATH_results, "BP_UP-barplot.pdf", sep=""), width = 6, height = 3)
mutate(ego2, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")
dev.off()

pdf(file = paste(PATH_results, "BP_UP-upset.pdf", sep=""), width = 15, height = 4)
upsetplot(ego2)
dev.off()

write.csv(ego2, paste(PATH_results, "BP_UP.csv"))


#--------

# Cellular Compartment
# NOTE: empty for UP
ego <- enrichGO(gene          = test_df,
                universe      = background,
                OrgDb         = org.Mm.eg.db,
                keyType       = "SYMBOL",
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

# remove redundancy in the GO terms
ego2    <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min, measure = 'Wang')

pdf(file = paste(PATH_results, "Compartment_DOWN-network.pdf", sep=""), width = 8, height = 8)
goplot(ego2)
dev.off()

pdf(file = paste(PATH_results, "Compartment_DOWN-barplot.pdf", sep=""), width = 6, height = 3)
mutate(ego2, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")
dev.off()

pdf(file = paste(PATH_results, "Compartment_DOWN-upset.pdf", sep=""), width = 15, height = 4)
upsetplot(ego2)
dev.off()

write.csv(ego2, paste(PATH_results, "Compartment_DOWN.csv"))

#--------

# Molecular Function
ego <- enrichGO(gene          = test_df,
                universe      = background,
                OrgDb         = org.Mm.eg.db,
                keyType       = "SYMBOL",
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

# remove redundancy in the GO terms
ego2    <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min, measure = 'Wang')

pdf(file = paste(PATH_results, "MF_UP-network.pdf", sep=""), width = 8, height = 8)
goplot(ego2)
dev.off()

pdf(file = paste(PATH_results, "MF_UP-barplot.pdf", sep=""), width = 6, height = 3)
mutate(ego2, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")
dev.off()

pdf(file = paste(PATH_results, "MF_UP-upset.pdf", sep=""), width = 15, height = 4)
upsetplot(ego2)
dev.off()

write.csv(ego2, paste(PATH_results, "MF_UP.csv"))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# check bias in receptor types for what is affected by chemo
results <- read.csv("./output/explants/ DEP-analysis-limma_drugeffect.csv", row.names = 1)
sig_de  <-  read.csv("./output/explants/ DEP-analysis-limma_drugeffect_sig.csv", row.names = 1)

# Check type of protein regulated
pg_types  <- readxl::read_excel("./data/interactome_list_v3.1_large.xlsx")
rec_types <- data.frame(symbol = pg_types$Rec_symbol, type = pg_types$Rec_type)
lig_types <- data.frame(symbol = pg_types$Lig_symbol, type = pg_types$Lig_type)

pg_types <- rbind(rec_types, lig_types)

head(pg_types)

# Convert interactome list to mouse identifiers using biomart
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",  host = "https://dec2021.archive.ensembl.org/")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",  host = "https://dec2021.archive.ensembl.org/")
  
  genesV2 <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x, mart = human,
                    attributesL = c("mgi_symbol"), martL = mouse, uniqueRows = TRUE)
  
  humanx <- unique(genesV2)
  return(humanx)
}

# biomaRt currently down (22.01.2025)
# genelist <- pg_types$symbol #extract symbols
# genelist <- convertHumanGeneList(genelist)
# 
# pg_types <- merge(genelist, pg_types, 
#                   by.x = "HGNC.symbol", 
#                   by.y = "symbol", 
#                   all.x = TRUE)

ortholgoues <- read.csv("../human_mouse_1to1_orthologs.csv")

pg_types <- merge(ortholgoues, pg_types, 
                  by.x = "human", 
                  by.y = "symbol", 
                  all.y = TRUE)


pg_types <- pg_types %>% distinct() # remove duplicates

head(pg_types)

write.csv(pg_types, "./data/published/receptortypes.csv")

# filter data of interest
sig_de <- sig_de[abs(sig_de$logFC) > 0.5, ] 

type_deps <- pg_types[pg_types$mouse %in% rownames(sig_de), ]
type_back <- pg_types[pg_types$mouse %in% rownames(results), ]

# type_deps <- pg_types[pg_types$MGI.symbol %in% rownames(sig_de), ]
# type_back <- pg_types[pg_types$MGI.symbol %in% rownames(results), ]

type_deps <- as.data.frame(table(type_deps$type))
type_back <- as.data.frame(table(type_back$type))

table(sig_de$adj.P.Val < 0.05)
table(results$adj.P.Val < 0.05)

type_deps <- type_deps %>% mutate(Proportion = Freq / 672)
type_back <- type_back %>% mutate(Proportion = Freq / 2928)

type_deps$Source <- "DEP"
type_back$Source <- "General"

combined_counts  <- rbind(type_deps, type_back)
head(combined_counts)

g <- ggplot(combined_counts, aes(x = Var1, y = Proportion, fill = Source)) 
g <- g + geom_bar(stat = "identity", position = "dodge") 
g <- g + theme_bw() +
  labs(title = "Receptor Bias from Oxaliplatin treatment",
       x = "",
       y = "Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
g <- g + scale_fill_manual(values = c("DEP" = "#a97cde", "General" = "#6ec083"),
    name = "Source")

print(g)

PATH_results = "./output/explants/"

pdf(paste0(PATH_results, "receptor_proportions.pdf"), height = 5, width = 12)
print(g)
dev.off()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# plot data by ion channels
df      <- read.csv("./data/explants/combined_genes.csv", row.names = 1, check.names = FALSE)
colData <- read.csv("./data/explants/colData.csv", row.names = 1)

enrichments  <- read.csv("./output/explants/tissue-enr/enrichments.csv", row.names = 1)

ion.channels <- read.csv("./data/published/receptortypes.csv", header = TRUE) #from hDRG prot paper

# Extract neuronal enriched genes (see 06-explant-enrichment.R)

rownames(df) <- trimws(rownames(df))
mat <- df[rownames(df) %in% enrichments$Gene, ]

head(mat)

colData <- colData[colData$sampleID %in% colnames(mat), ]
index   <- match(colnames(mat), colData$sampleID)
colData <- colData[index, ]

head(colnames(mat))
head(colData$sampleID)

ion.channels <- na.omit(ion.channels)
ion.channels <- ion.channels$MGI.symbol[ion.channels$type == "ion channel"]

ion.channels <- trimws(ion.channels)

rownames(mat) <- trimws(as.character(rownames(mat)))
data <- mat[rownames(mat) %in% ion.channels, ]

scaled_expression <- t(scale(t(data), center = TRUE))

match_index <- match(colnames(scaled_expression), colData$sampleID)
colData_reordered <- colData[match_index, ]

tissue_list <- as.factor(paste(colData_reordered$Condition,"-",colData_reordered$Turbo))

# Create a color mapping for colData

tissue_colors <- c("Ox - T" = "#3b92df",
                   "Ox - TC" = "#bed1e1",
                   "V - T" = "#edb127",
                   "V - TC" = "#f7e1ae"
)

# Assign colors to colData levels
col_fun <- tissue_colors[tissue_list]

scaled_expression <- t(scale(t(data)))
scaled_expression[is.na(scaled_expression)] <- 0

#----------

# Make a fresh colour gradiant so missing values stand out
min_val <- min(scaled_expression, na.rm = TRUE)
max_val <- max(scaled_expression, na.rm = TRUE)
viridis_colors <- viridis(100)

col_fun2 <- colorRamp2(
  c(min_val, -0.0000001, 0, 0.0000001, max_val),  # Data range with zero explicitly included
  c(viridis_colors[50], "white", "grey", "white", viridis_colors[1]))

#---------

# Create Heatmap
ht_list <- ComplexHeatmap::Heatmap(scaled_expression,
                                   #name = "Expression",
                                   col= col_fun2,
                                   clustering_distance_columns = "manhattan",
                                   #cluster_rows = TRUE,
                                   #cluster_columns = TRUE,
                                   show_row_names = TRUE,
                                   show_column_names = FALSE, #set to TRUE to double check colour legend
                                   row_title = "Proteins",
                                   row_dend_side = "left",
                                   top_annotation = HeatmapAnnotation(tissue = tissue_list, col = list(tissue = col_fun))
)

draw(ht_list, heatmap_legend_side = "right")

PATH_results = "./output/"

pdf(file = paste0(PATH_results, "ionchannel-heatmap.pdf"), height = 5, width = 6)
draw(ht_list, heatmap_legend_side = "right")
dev.off()

#-------------------------------------------------------------------------------

# TurboID explants vs WCL (see 06-explants-wcl.R)
dir.create("./output/explants/turboVwcl")
PATH_results = "./output/explants/turboVwcl/"

# df      <- read.csv("./data/explants/combined_genes.csv", row.names = 1, check.names = FALSE)
# colData <- read.csv("./data/explants/colData.csv", row.names = 1)
# enrichments  <- read.csv("./output/explants/tissue-enr/enrichments.csv", row.names = 1)

results.turbo <- read.csv("./output/explants/ DEP-analysis-limma_drugeffect.csv", row.names = 1)
results.wcl   <- read.csv("./output/explants/wcl/ DEP-analysis-limma_ox.csv", row.names = 1)

turbo <- results.turbo[results.turbo$adj.P.Val < 0.05, ]
wcl   <- results.wcl[results.wcl$adj.P.Val < 0.05, ]

merged_data <- merge(turbo, wcl, by = "row.names", all = FALSE, suffixes = c(".turbo",".wcl"))
head(merged_data)

correlation <- cor(merged_data$logFC.turbo, merged_data$logFC.wcl, use = "complete.obs")

g <- ggplot(merged_data, aes(x = logFC.turbo, y = logFC.wcl)) 
g <- g + geom_point(color = "#9357f1") 
g <- g + geom_smooth(method = "lm", color = "grey", se = TRUE) 
g <- g + labs(title = paste("DEPs, Turbo vs WCL
Corr=",round(correlation, 2)),
              x = "logFC, Turbo",
              y = "logFC, WCL") +
  theme_bw()

print(g)

pdf(file = paste0(PATH_results, "correlation_LFCp0.05.pdf"), height = 4, width = 4)
print(g)
dev.off()

write.csv(merged_data, "./output/explants/turboVwcl/DEPs_across_datasets.csv")

#-----------------------

turbo <- results.turbo
wcl   <- results.wcl

merged_data <- merge(turbo, wcl, by = "row.names", all = FALSE, suffixes = c(".turbo",".wcl"))
head(merged_data)

correlation <- cor(merged_data$AveExpr.turbo, merged_data$AveExpr.wcl, use = "complete.obs")

g <- ggplot(merged_data, aes(x = AveExpr.turbo, y = AveExpr.wcl)) 
g <- g + geom_point(color = "#9357f1") 
g <- g + geom_smooth(method = "lm", color = "grey", se = TRUE) 
g <- g + labs(title = paste("AveExpr (DEPs), Turbo vs WCL
Corr=",round(correlation, 2)),
              x = "Expression, Turbo",
              y = "Expression, WCL") +
  theme_bw()

print(g)

pdf(file = paste0(PATH_results, "correlation_AvgExp_p0.05.pdf"), height = 4, width = 4)
print(g)
dev.off()

#-------------------------------------------------------------------------------

overlap <- read.csv("./output/explants/turboVwcl/DEPs_across_datasets.csv", row.names = 1)
overlap <- na.omit(overlap)
rownames(overlap) <- overlap$Row.names
overlap$Row.names <- NULL

head(overlap)

overlap$logFC_overlap <- ifelse(
  sign(overlap$logFC.turbo) == sign(overlap$logFC.wcl) & overlap$logFC.turbo > 0, "+/+",
  ifelse(sign(overlap$logFC.turbo) == sign(overlap$logFC.wcl) & overlap$logFC.turbo < 0, "-/-",
    "Mismatched"))

head(overlap)

df      <- read.csv("./data/explants/combined_genes.csv", row.names = 1, check.names = FALSE)
colData <- read.csv("./data/explants/colData.csv", row.names = 1)
enrichments  <- read.csv("./output/explants/tissue-enr/enrichments.csv", row.names = 1)

mat <- df[rownames(df) %in% enrichments$Gene, ]
head(colData$sampleID)

data <- mat[rownames(mat) %in% rownames(overlap), ]
# data <- na.omit(data)

scaled_expression <- t(scale(t(data), center = TRUE))

match_index <- match(colnames(scaled_expression), colData$sampleID)
colData_reordered <- colData[match_index, ]

tissue_list <- as.factor(paste(colData_reordered$Condition,"-",colData_reordered$Turbo))

# Create a color mapping for colData

tissue_colors <- c("Ox - T" = "#3b92df",
                   "Ox - TC" = "#bed1e1",
                   "V - T" = "#edb127",
                   "V - TC" = "#f7e1ae"
)

# Assign colors to colData levels
col_fun <- tissue_colors[tissue_list]

scaled_expression <- t(scale(t(data)))
scaled_expression[is.na(scaled_expression)] <- 0

# Make a fresh colour gradiant so missing values stand out
min_val <- min(scaled_expression, na.rm = TRUE)
max_val <- max(scaled_expression, na.rm = TRUE)
viridis_colors <- viridis(100)

col_fun2 <- colorRamp2(
  c(min_val, -0.0000001, 0, 0.0000001, max_val),  # Data range with zero explicitly included
  c(viridis_colors[50], "white", "grey", "white", viridis_colors[1]))

#---------

# Create Heatmap
ht_list <- ComplexHeatmap::Heatmap(scaled_expression,
                                   #name = "Expression",
                                   col= col_fun2,
                                   clustering_distance_columns = "manhattan",
                                   cluster_rows = TRUE,
                                   cluster_columns = TRUE,
                                   show_row_names = FALSE,
                                   show_column_names = FALSE, #set to TRUE to double check colour legend
                                   row_title = "Proteins",
                                   row_dend_side = "left",
                                   top_annotation = HeatmapAnnotation(tissue = tissue_list, col = list(tissue = col_fun)))


draw(ht_list, heatmap_legend_side = "right")

pdf(file = paste0(PATH_results, "matched_DEP_heatmap.pdf"), height = 5, width = 6)
draw(ht_list, heatmap_legend_side = "right")
dev.off()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Comparison to other chemo DRG data?
# Proteomic: https://doi.org/10.1016/j.jprot.2022.104682
# translatome: https://doi.org/10.1523/JNEUROSCI.2661-18.2018

results     <- read.csv("./output/explants/ DEP-analysis-limma_drugeffect.csv", row.names = 1)
translatome <- read.csv("./data/megat2018_paclicatalDEGs.csv") #combined from Fig4-2 and 4-3 from Megat 2018, using their definition of DEG

results <- results[results$adj.P.Val < 0.05, ]

head(translatome)
head(results)

common_genes <- intersect(rownames(results), translatome$Gene.name)

# Step 2: Subset logFC values for matched genes
results_matched <- results[common_genes, "logFC", drop=FALSE]
translatome_matched <- translatome[translatome$Gene.name %in% common_genes, c("Gene.name", "Log2FC")]

# Step 3: Merge datasets
merged_data <- merge(
  data.frame(Gene.name = rownames(results_matched), logFC = results_matched$logFC),
  translatome_matched,
  by = "Gene.name"
)

merged_data <- na.omit(merged_data)

# Step 4: Calculate correlation
correlation <- cor(merged_data$logFC, merged_data$Log2FC)

# Step 5: Plot
g <- ggplot(merged_data, aes(x = logFC, y = Log2FC)) 
g <- g + geom_point(color = "#9357f1") 
g <- g + geom_smooth(method = "lm", color = "grey", se = TRUE) 
g <- g + labs(title = paste("DEPs, vs Megat 2018
Corr=",round(correlation, 2)),
x = "logFC (TurboID)",
y = "logFC (Megat2018)"
) +
  theme_bw()

print(g)

pdf(file = paste0(PATH_results, "Megat2018_correlation_turbo.pdf"), height = 4, width = 4)
print(g)
dev.off()

#-------------------------

yang <- read.csv("./data/published/Yang2022_OxaVsVeh_shortterm_stable2.csv")

head(yang)
yang <- yang[, colnames(yang) %in% c("Gene_Symbol", "log2FC", "log2.protein.intensity", "Adjusted_P_value"), ]

head(results)

results.turbo <- results

# # select only DEPs from each
inhouse <- results.turbo[results.turbo$adj.P.Val < 0.05, ]
yang    <- yang[yang$Adjusted_P_value < 0.05, ]

merged_data <- merge(inhouse, yang, by.x = "row.names", by.y = "Gene_Symbol", all = FALSE)
head(merged_data)

#intuitive naming for saved csv
colnames(merged_data) <- c("GeneID", "explant.logFC", "explant.AveExpr", "explant.t", 
                           "explant.P.Value", "explant.adj.P.Val", "explant.B",
                           "yang.logFC", "yang.AveExpr", "yang.adj.P.Val")

correlation <- cor(merged_data$explant.logFC, merged_data$yang.logFC, use = "complete.obs")

g <- ggplot(merged_data, aes(x = explant.logFC, y = yang.logFC)) 
g <- g + geom_point(color = "#9357f1") 
g <- g + geom_smooth(method = "lm", color = "grey", se = TRUE) 
g <- g + labs(title = paste("DEPs, vs Yang et al 2022.
Corr=",round(correlation, 2)),
              x = "logFC (TurboID)",
              y = "logFC (Yang2022)"
) +
  theme_bw()

print(g)

pdf(file = paste0(PATH_results, "yang2022_correlation_p0.05_turbo.pdf"), height = 4, width = 4)
print(g)
dev.off()

