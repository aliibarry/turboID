library(readxl)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggbiplot)
library(ggrepel)
library(viridis)
library(circlize)
library(ComplexHeatmap)
library(matrixStats)
library(gridExtra)
library(stringr)

dir.create("./output/wcl")
PATH_results = "./output/wcl/"

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

df <- read.csv("./data/20240614_wcl_combined.csv", check.names = FALSE)

column_names <- colnames(df)[-1]
metadata <- do.call(rbind, strsplit(column_names, "_"))

metadata <- as.data.frame(metadata)
colnames(metadata) <- c("Sample", "Tissue", "Type", "date") #clarify TC/T and ID from Julia
metadata$sampleID <- column_names

#-------------------------------------------------------------------------------

rownames(df) <- df$GeneID
df$GeneID <- NULL

colData <- metadata

df <- df[, colnames(df) %in% colData$sampleID[colData$Type != "pool" | colData$Tissue == "paw"]]
colData <- colData[colData$Type != "pool" | colData$Tissue == "paw", ]
colData$Tissue[colData$Tissue == "LDRG"] <- "DRG"
colData$Tissue[colData$Tissue == "dLSC"] <- "LSC"

# Reshape the dataframe to long format
df_long <- df %>%
  pivot_longer(cols = everything(), names_to = "sampleID", values_to = "Expression") %>%
  left_join(colData, by = "sampleID")

summed_expression <- df_long %>%
  group_by(sampleID, Tissue) %>%
  summarize(SumCount = sum(Expression, na.rm = TRUE), .groups = "drop") %>%
  mutate(Tissue = factor(Tissue, levels = c("LSC", "DRG", "SCN", "paw")))

mean_intensity <- df_long %>%
  group_by(sampleID, Tissue) %>%
  summarize(MeanIntensity = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
  mutate(Tissue = factor(Tissue, levels = c("LSC", "DRG", "SCN", "paw")))

detected_counts <- df_long %>%
  group_by(sampleID, Tissue) %>%
  summarize(DetectedCount = sum(!is.na(Expression)), .groups = "drop") %>%
  mutate(Tissue = factor(Tissue, levels = c("LSC", "DRG", "SCN", "paw")))

#-------

g1 <- ggplot(summed_expression, aes(x = Tissue, y = SumCount, fill = Tissue)) 
g1 <- g1 + geom_boxplot() 
g1 <- g1 + theme_classic() +
  labs(title = "", 
       x = "", 
       y = "Summed Expression")

print(g1)

#------

g2 <- ggplot(mean_intensity, aes(x = Tissue, y = MeanIntensity, fill = Tissue)) 
g2 <- g2 + geom_boxplot() 
g2 <- g2 + theme_classic() +
  labs(title = "", 
       x = "", 
       y = "Mean Intensity")

print(g2)

#-------

g3 <- ggplot(detected_counts, aes(x = Tissue, y = DetectedCount, fill = Tissue)) 
g3 <- g3 + geom_boxplot() 
g3 <- g3 + theme_classic() +
  labs(title = "", 
       x = "", 
       y = "PG Count")

print(g3)


#------

pdf(file = paste(PATH_results, "expression_sum.pdf", sep="/"), width = 5, height = 4)
print(g1)
dev.off()

pdf(file = paste(PATH_results, "expression_mean.pdf", sep="/"), width = 5, height = 4)
print(g2)
dev.off()

pdf(file = paste(PATH_results, "PG_counts.pdf", sep="/"), width = 5, height = 4)
print(g3)
dev.off()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Heatmap of neuronal/myelin genes on interest
neurons <- read.csv("./data/neuronal-genes_mouse.csv", header = FALSE) #from hDRG prot paper

# data <- data %>%
#   separate(genes, into = paste0("genes"), sep = ";", remove = TRUE) 

df$GeneID <- rownames(df)
df$GeneID <- trimws(as.character(df$GeneID))

data <- df %>%
  mutate(GeneID = gsub("\\n", "", GeneID)) %>%
  distinct(GeneID, .keep_all = TRUE)

# data$genes <- NULL
# expression_data <- data[complete.cases(data),]

data <- data[data$GeneID %in% neurons$V1, ]

head(data)
rownames(data) <- data$GeneID
data$GeneID <- NULL

data <- as.matrix(data)

scaled_expression <- t(scale(t(data), center = TRUE))

metadata <- colData

match_index <- match(colnames(scaled_expression), metadata$sampleID)
metadata_reordered <- metadata[match_index, ]

tissue_list <- as.factor(metadata_reordered$Tissue)

# Create a color mapping for metadata
tissue_colors <- c("DRG" = "#f55c3a",
                   "LSC" = "#edb127",
                   "SCN" = "#29c99d",
                   "paw" = "#95459b"
                   )

# Assign colors to metadata levels
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
                   cluster_rows = TRUE,
                   cluster_columns = TRUE,
                   show_row_names = TRUE,
                   show_column_names = FALSE, #set to TRUE to double check colour legend
                   row_title = "Proteins",
                   row_dend_side = "left",
                   top_annotation = HeatmapAnnotation(tissue = tissue_list, col = list(tissue = col_fun))
)

draw(ht_list, heatmap_legend_side = "right")

PATH_results = "./output/wcl/"

pdf(file = paste0(PATH_results, "/neuronal-heatmap.pdf"), height = 7, width = 6)
draw(ht_list, heatmap_legend_side = "right")
dev.off()

df[sapply(df, is.numeric)] <- log2(df[sapply(df, is.numeric)])

write.csv(df, "./data/wcl-matrix.csv")
write.csv(colData, "./data/wcl-colData.csv")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

df      <- read.csv("./data/wcl-matrix.csv", header = TRUE)
colData <- read.csv("./data/wcl-colData.csv", header = TRUE)

head(df)

mat <- as.matrix(df[, which(colnames(df) %in% colData$sampleID)])
mat <- mat[complete.cases(mat), ]

rv     <- matrixStats::rowVars(mat) # calculate variance per row (ie. per gene)
select <- order(rv, decreasing=TRUE)[seq_len(min(2000, length(rv)))]

pca <- prcomp(t(mat[select,]), center = TRUE, scale. = TRUE)
pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)

index   <- match(colnames(mat), colData$sampleID)
colData <- colData[index, ]

tissue <- as.factor(colData$Tissue)

g1 <- ggbiplot(pca, choices = c(1,2),
               groups = interaction(tissue),
               #ellipse = TRUE,
               ellipse.prob = 0.95,
               labels = NULL,
               point.size = 4,
               labels.size = 4, alpha = 1, var.axes = FALSE,
               circle  = TRUE, circle.prob = 0.5,
               varname.size = 3,
               varname.adjust = 1.5,
               varname.abbrev = FALSE)

g1 <- g1 + theme(legend.position = 'right', aspect.ratio= 1) + theme_bw()

print(g1)

PATH_results = "./output/wcl/"

pdf(file = paste(PATH_results, "pca.pdf", sep=""), width = 4, height = 4)
print(g1)
dev.off()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Merge Turbo + WCL together for analyses
wcl_df   <- read.csv("./data/wcl-matrix.csv", header = TRUE, row.names = 1)
wcl_meta <- read.csv("./data/wcl-colData.csv", header = TRUE, row.names = 1)

turbo_df   <- read.csv("./data/matrix-for-limma.csv", header = TRUE)
turbo_meta <- read.csv("./data/colData-for-limma.csv", header = TRUE)

to_keep <- c("Tissue", "Turbo", "sampleID", "Dataset")
wcl_meta$Turbo <- "NA"

wcl_meta$Dataset   <- "WCL"
turbo_meta$Dataset <- "Turbo"

wcl_meta   <- wcl_meta[, colnames(wcl_meta) %in% to_keep]
turbo_meta <- turbo_meta[, colnames(turbo_meta) %in% to_keep]

head(wcl_meta)
head(turbo_meta)

colData <- rbind(wcl_meta, turbo_meta)

wcl_df   <- as.data.frame(wcl_df) %>%
  mutate(GeneID = sub(";.*", "", GeneID)) %>%
  distinct(GeneID, .keep_all = TRUE)

turbo_df <- as.data.frame(turbo_df) %>%
  distinct(genes, .keep_all = TRUE)

turbo_df$GeneID <- turbo_df$genes
turbo_df$proteins <- NULL
turbo_df$genes    <- NULL

# numeric_cols <- setdiff(names(turbo_df), "GeneID")
# turbo_df[, numeric_cols] <- scale(turbo_df[, numeric_cols], center = TRUE, scale = TRUE)
# 
# numeric_cols <- setdiff(names(wcl_df), "GeneID")
# wcl_df[, numeric_cols] <- scale(wcl_df[, numeric_cols], center = TRUE, scale = TRUE)

df <- merge(wcl_df, turbo_df, by = "GeneID", all = TRUE)

rownames(df) <- df$GeneID
df$GeneID <- NULL

#-------------------------------------------------------------------------------

# Look at the enrichment of neuronal genes in the paw
head(colData)

sub_meta <- colData[colData$Tissue == "paw",]
sub_df   <- df[, colnames(df) %in% sub_meta$sampleID]

neurons <- read.csv("./data/neuronal-genes_mouse.csv", header = FALSE) #from hDRG prot paper

sub_df$GeneID <- rownames(sub_df)
sub_df$GeneID <- trimws(as.character(sub_df$GeneID))

data <- sub_df %>%
  mutate(GeneID = gsub("\\n", "", GeneID)) %>%
  distinct(GeneID, .keep_all = TRUE)

data <- data[data$GeneID %in% neurons$V1, ]

head(data)
rownames(data) <- data$GeneID
data$GeneID <- NULL

data <- as.matrix(data)
data <- data[rowSums(is.na(data)) < ncol(data), ]

scaled_expression <- data
scaled_expression[is.na(scaled_expression)] <- 0

match_index <- match(colnames(scaled_expression), sub_meta$sampleID)
metadata_reordered <- sub_meta[match_index, ]

tissue_list <- as.factor(paste(metadata_reordered$Tissue,metadata_reordered$Dataset, sep="."))

# Create a color mapping for metadata
tissue_colors <- c("paw.WCL" = "#dab0de",
                   "paw.Turbo" = "#95459b"
)

# Assign colors to metadata levels
col_fun <- tissue_colors[tissue_list]

#----------

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

# for unscaled data (log2)
# min_val <- min(scaled_expression, na.rm = TRUE)
# max_val <- max(scaled_expression, na.rm = TRUE)
# viridis_colors <- viridis(100)
# 
# col_fun2 <- colorRamp2(
#   c(0, min_val, max_val),  # Data range with zero explicitly included
#   c("grey", "white", viridis_colors[50]))

#---------

# Create Heatmap
ht_list <- ComplexHeatmap::Heatmap(scaled_expression,
                                   #name = "Expression",
                                   col= col_fun2,
                                   clustering_distance_columns = "manhattan",
                                   cluster_rows = TRUE,
                                   cluster_columns = TRUE,
                                   show_row_names = TRUE,
                                   show_column_names = FALSE, #set to TRUE to double check colour legend
                                   row_title = "Proteins",
                                   row_dend_side = "left",
                                   top_annotation = HeatmapAnnotation(tissue = tissue_list, col = list(tissue = col_fun))
)

draw(ht_list, heatmap_legend_side = "right")

PATH_results = "./output/wcl/"

pdf(file = paste0(PATH_results, "/neuronal-heatmap_paw.pdf"), height = 5, width = 6)
draw(ht_list, heatmap_legend_side = "right")
dev.off()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(limma)
library(readxl)
library(ggplot2)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)

PATH_results = "./output/enrichments/subtype/"

wcl_df      <- read.csv("./data/wcl-matrix.csv", header = TRUE, row.names = 1)
wcl_meta    <- read.csv("./data/wcl-colData.csv", header = TRUE, row.names = 1)
enrichments <- read.csv("./output/enrichments_75filt.csv", check.names = FALSE, header = TRUE, row.names = 1)

zheng_full  <- read.csv("./data/zheng_gl.csv", row.names = 1)
zheng_full  <- zheng_full[, colnames(zheng_full) %in% c("gs", "symbol")]

#--------------------------

# Build background of WCL within tissue + enriched proteins for that tissue
toi = "SCN" #tissue of interest

samples <- wcl_meta$sampleID[wcl_meta$Tissue == toi]
background_df <- wcl_df[, colnames(wcl_df) %in% samples]

# remove columns within any valid values
index <- apply(background_df, 1, function(x) all(is.na(x)))
background_df <- background_df[ !index, ]

# WCL background
background_gl1 <- data.frame(gs = "background",
                            symbol = rownames(background_df))

# Turbo background
background_gl2 <- data.frame(gs = "background",
                             symbol = enrichments$Gene[enrichments$Tissue == toi])

background_gl <- rbind(background_gl1, background_gl2, keep.all = TRUE) ## combined background

# merge background with rest of gene sets, to bypass bug to assign universe
zheng_gl <- zheng_full[zheng_full$symbol %in% background_gl$symbol, ]
gl_full  <- rbind(zheng_gl, background_gl, keep.all = TRUE)

#---------                               

# test overenrichment of subpopulations
test_df  <- enrichments$Gene[enrichments$Tissue == toi] 

options(enrichment_force_universe = TRUE)
ego <- enricher(gene          = test_df,
                # universe      = background$Gene, #doesn't work, added bg in as a fake geneset to include it
                minGSSize     = 5,
                maxGSSize     = 2000,
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                pAdjustMethod = "BH",
                TERM2GENE     = gl_full)

result_df <- data.frame(ego@result)

result_df <- result_df %>%
  separate(GeneRatio, into = c("Gene_in_Set", "Total_Genes"), sep = "/", convert = TRUE) %>%
  separate(BgRatio, into = c("Bg_in_Set", "Total_Bg"), sep = "/", convert = TRUE) %>%
  mutate(GeneRatio_numeric = Gene_in_Set / Bg_in_Set,
         BgRatio_numeric = Bg_in_Set / Total_Bg)

# Print results
head(result_df)

write.csv(result_df, file = paste0(PATH_results, "/ORA_Zheng_",toi,"-2.csv"))

g <- ggplot(result_df, aes(x=(Description), y=GeneRatio_numeric, colour=p.adjust, size=Count))
g <- g + geom_point() + theme_bw() + ggtitle(toi) +
  theme(axis.text.y = element_text(size= 12, colour= "black", hjust = 1), 
        axis.text.x = element_text(size=10, angle = 45, hjust= 1), 
        legend.text = element_text(size=10), 
        axis.title.x = element_blank(),
        plot.title=element_text(size=rel(1), hjust = 1)) +
  theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm")) +  
  labs(y="GS Proportion", colour="p.adj", size="Count") +
  scale_colour_gradient(limits = c(0, 1)) +
  scale_size_continuous(limits = c(10, 500)) +
  ylim(0, 0.75) 

pdf(paste0(PATH_results, "ORA_Zheng_",toi,".pdf"), height = 4, width = 4)
print(g)
dev.off()

print(g)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------     
#-------------------------------------------------------------------------------

# SCRATCH, LFC cann't account for sample prep differences
# enrichments <- read.csv("./output/enrichments_75filt.csv", 
#                         check.names = FALSE, header = TRUE, row.names = 1)
# head(enrichments)
# 
# df_sub <- df[rownames(df) %in% enrichments$Gene, ]
# 
# head(df_sub)
# 
# colData <- colData[colData$Turbo  != "TC", ]
# colData$Tissue.Set <- as.factor(paste(colData$Tissue,colData$Dataset, sep = "."))
# 
# #remove turbo control samples
# df_sub <- df_sub[, colnames(df_sub) %in% colData$sampleID]
# 
# mat <- df_sub
# mat <- mat[,colnames(mat) %in% colData$sampleID]
# 
# index   <- match(colnames(mat), colData$sampleID)
# colData <- colData[index, ]
# 
# head(colData)
# mat[1:5, 1:5]
# 
# design <- model.matrix(~ 0 + Tissue.Set, data = colData) #select design
# head(design)
# 
# fit <- lmFit(mat, design)
# 
# colnames(fit) #check possible comparisons
# 
# # Run hypothesis testing, adjust tissue as needed. Only for ranked LFc. Data not normalized between WCL + Turbo
# contrast_matrix <- makeContrasts( Tissue.SetSCN.Turbo - Tissue.SetSCN.WCL, levels = design)
# fit2 <- contrasts.fit(fit, contrast_matrix)
# fit2 <- eBayes(fit2)
# 
# results <- topTable(fit2, adjust="BH", number=Inf)
# 
# head(results)
# 
# sig_de <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, , drop = FALSE] #extract significant proteins if there
# head(sig_de)
# 
# sig_de <- na.omit(sig_de)
# dim(sig_de)
# 
# # volcano plotting
# res       <- as.data.frame(results) 
# mutateddf <- mutate(res, Sig=ifelse(res$adj.P.Val<0.05 & abs(res$logFC)>1, "FDR < 0.05 & LFC > 1", ifelse("NS")))
# input     <- cbind(gene=rownames(res), mutateddf) 
# 
# volc = ggplot(input, aes(logFC, -log10(P.Value))) + geom_point(aes(col=Sig)) +
#   scale_color_manual(values = c("#B63679ff", "grey")) + 
#   ggrepel::geom_text_repel(data=subset(input, input$gene %in% rownames(sig_de)),
#                            aes(label=gene), size=4, segment.alpha= 0.2, force =2, max.overlaps=16) 
# volc <- volc + theme_bw() + theme(aspect.ratio=1)
# volc <- volc + theme(legend.position="bottom", axis.text.y = element_text(size= 12, face="bold"), 
#                      axis.title.y = element_text(size=14), axis.title.x = element_text(size= 14), 
#                      axis.text.x = element_text(size= 12), legend.title=element_text(size=14), 
#                      legend.text=element_text(size=14), plot.title=element_text(size=12, hjust = 0.5)) + 
#   ggtitle("SCN")
# print(volc)
# 
# 
# table(sig_de$logFC > 1)
# 
# write.csv(results, paste(PATH_results, "TurboLFC-SCN.csv"))
# 
# #-------------------------------------------------------------------------------
# 
# load("../bulkseq/DE/zheng2019_genelists.RData")
# ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# gene_map <- getBM(
#   attributes = c("ensembl_gene_id", "mgi_symbol"),
#   filters = "ensembl_gene_id",
#   values = rownames(zheng_gl),
#   mart = ensembl
# )
# gene_map <- gene_map[!duplicated(gene_map$ensembl_gene_id), ]
# rownames(gene_map) <- gene_map$ensembl_gene_id
# zheng_gl$symbol <- gene_map[rownames(zheng_gl), "mgi_symbol"]
# 
# write.csv(zheng_gl, "./data/zheng_gl.csv")
# 
# # Load tissue of interest
# results <- read.csv("./output/wcl/ TurboLFC-LSC.csv", row.names = 1)
# 
# lfc_vector <- results$logFC
# names(lfc_vector) <- rownames(results)
# 
# lfc_vector <- sort(lfc_vector, decreasing = TRUE) # Sort log2 fold change values in descending order
# 
# set.seed(52)
# 
# gsea_results <- GSEA(
#   geneList = lfc_vector,
#   minGSSize = 0,
#   maxGSSize = 2000,
#   eps = 0,
#   nPermSimple = 10000,
#   pvalueCutoff = 1, 
#   seed = TRUE, 
#   pAdjustMethod = "BH",
#   TERM2GENE = dplyr::select(zheng_gl, gs, symbol) 
# )
# 
# #head(gsea_results@result)
# gsea_result_df <- data.frame(gsea_results@result)
# write.csv(gsea_result_df, file = paste0(PATH_results, "/GSEA_Zheng_LSC.csv"))
# 
# g <- ggplot(gsea_result_df, aes(x=(Description), y=NES, colour=p.adjust, size=setSize))
# g <- g + geom_point() + theme_bw() + ggtitle("LSC") +
#   theme(axis.text.y = element_text(size= 12, colour= "black", hjust = 1), 
#         axis.text.x = element_text(size=10, angle = 45, hjust= 1), 
#         legend.text = element_text(size=10), 
#         axis.title.x = element_blank(),
#         plot.title=element_text(size=rel(1), hjust = 1)) +
#   theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm")) +  
#   labs(y="Enrichment Score", colour="p value", size="Count")
# 
# pdf(paste0(PATH_results, "GSEA_LSC.pdf"), height = 4, width = 4)
# print(g)
# dev.off()