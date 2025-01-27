library(dplyr)
library(tidyr)
library(diann)
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

dir.create("./output/explants/wcl")
PATH_results = "./output/explants/wcl/"

# update to report without pools once run
df <- diann::diann_load("./data/explants/wcl/DRGexpl_report.tsv")

precursors <- diann_matrix(df, q = 0.01)

peptides.maxlfq <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,], 
                                sample.header = "Run",
                                group.header="Stripped.Sequence", 
                                id.header = "Precursor.Id", 
                                quantity.header = "Precursor.Normalised")

gene.groups <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$GG.Q.Value <= 0.01,], 
                            sample.header = "Run",
                            group.header="Genes", 
                            id.header = "Precursor.Id", 
                            quantity.header = "Precursor.Normalised")

write.csv(precursors,      "./data/explants/wcl/precursors.csv", na="NA",eol = "\n", row.names = T)
write.csv(peptides.maxlfq, "./data/explants/wcl/peptides.csv", na="NA", eol = "\n", row.names = T)
write.csv(gene.groups,     "./data/explants/wcl/genegroups.csv", na="NA", eol = "\n", row.names = T)

#-------------------------------------------------------------------------------

mat <- peptides.maxlfq[complete.cases(peptides.maxlfq), ]
correlation_matrix <- cor(mat)

melted_corr_matrix <- melt(correlation_matrix)

# Create the ggplot
ComplexHeatmap::pheatmap(correlation_matrix,
                         main = "Correlation matrix",
                         border_color = "black",
                         cluster_rows = TRUE,
                         cluster_cols = TRUE,
                         col=viridis(100),
                         fontsize_row = 8,
                         fontsize_col = 8
)


pdf(file = paste(PATH_results, "correlation.pdf", sep=""), width = 5, height = 5)
ComplexHeatmap::pheatmap(correlation_matrix,
                         main = "Correlation matrix",
                         border_color = "black",
                         cluster_rows = TRUE,
                         cluster_cols = TRUE,
                         col=viridis(100),
                         fontsize_row = 8,
                         fontsize_col = 8
)

dev.off()

pre_counts  <- colSums(!is.na(precursors))
pre_counts  <- as.data.frame(pre_counts)

sorted_indices <- order(-pre_counts$pre_counts)
sorted_counts <- pre_counts$pre_counts[sorted_indices]
sorted_names <- rownames(pre_counts)[sorted_indices]

barplot(sorted_counts, 
        names.arg = sorted_names,
        #xlab = "Sample", 
        ylab = "Precursor Count", 
        main = "Counts by Sample",
        las = 2,  # Rotate x-axis labels by 45 degrees
        cex.names = 0.8)  # Decrease x-axis label text size


pdf(file = paste(PATH_results, "precursor-counts.pdf", sep=""), width = 8, height = 8)
par(mfrow=c(2,1))
barplot(sorted_counts, 
        names.arg = sorted_names,
        #xlab = "Sample", 
        ylab = "Precursor Count", 
        #main = "Counts by Sample",
        las = 2,
        cex.names = 0.8)
#abline(h = 15000, col = "red", lwd = 2)
dev.off()

# PG counts
pg_counts  <- colSums(!is.na(gene.groups))
pg_counts <- as.data.frame(pg_counts)

sorted_indices <- order(-pg_counts$pg_counts)
sorted_counts <- pg_counts$pg_counts[sorted_indices]
sorted_names <- rownames(pg_counts)[sorted_indices]

pdf(file = paste(PATH_results, "pg-counts.pdf", sep=""), width = 8, height = 8)
par(mfrow=c(2,1))
barplot(sorted_counts, 
        names.arg = sorted_names,
        #xlab = "Sample", 
        ylab = "PG Count", 
        #main = "Counts by Sample",
        las = 2,
        cex.names = 0.8)
#abline(h = 15000, col = "red", lwd = 2)
dev.off()

pdf(file = paste(PATH_results, "pg-boxplot.pdf", sep=""), width = 4, height = 4)
boxplot(sorted_counts)
dev.off()

#-------------------------------------------------------------------------------


df <- log2(gene.groups)

# Build corresponding metadata
sample_names  <- colnames(df)
sample_info   <- strsplit(sample_names, "_")
colData       <- as.data.frame(do.call(rbind, sample_info))
colnames(colData) <- c("Condition", "Concentration", "Slot", "Value", "Sample")
colData$Condition <- gsub("[0-9]", "", as.character(colData$Condition))

colData <- as.data.frame(cbind(sampleID = sample_names, colData))
head(colData)

# remove pooled samples
experimental <- colData$sampleID[colData$Condition == "WT"]
df <- df[,colnames(df) %in% experimental]
colData <- colData[colData$Condition == "WT", ]

# confirm matching order
df[1:5,1:5]
head(colData)

dim(df)

data <- as.data.frame(df) %>%
  mutate(GeneID = sub(";.*", "", rownames(df))) %>%  # Keep only the part before the semicolon
  distinct(GeneID, .keep_all = TRUE)

dim(data)

write.csv(data, "./data/explants/wcl/matrix.csv")
write.csv(colData, "./data/explants/colData.csv")

#-------------------------------------------------------------------------------

mat <- df[complete.cases(df), ]
correlation_matrix <- cor(mat)

melted_corr_matrix <- melt(correlation_matrix)

# Create the ggplot
ComplexHeatmap::pheatmap(correlation_matrix,
                         main = "Correlation matrix",
                         border_color = "black",
                         cluster_rows = TRUE,
                         cluster_cols = TRUE,
                         col=viridis(100),
                         fontsize_row = 8,
                         fontsize_col = 8
)

pdf(file = paste(PATH_results, "correlation.pdf", sep=""), width = 5, height = 5)
ComplexHeatmap::pheatmap(correlation_matrix,
                         main = "Correlation matrix, PGs",
                         border_color = "black",
                         cluster_rows = TRUE,
                         cluster_cols = TRUE,
                         col=viridis(100),
                         fontsize_row = 8,
                         fontsize_col = 8)
dev.off()

# PG counts
pg_counts  <- colSums(!is.na(df))
pg_counts  <- as.data.frame(pg_counts)

sorted_indices <- order(-pg_counts$pg_counts)
sorted_counts <- pg_counts$pg_counts[sorted_indices]
sorted_names <- rownames(pg_counts)[sorted_indices]

pdf(file = paste(PATH_results, "pg-counts.pdf", sep=""), width = 8, height = 8)
par(mfrow=c(2,1))
barplot(sorted_counts, 
        names.arg = sorted_names,
        #xlab = "Sample", 
        ylab = "PG Count", 
        #main = "Counts by Sample",
        las = 2,
        cex.names = 0.8)
#abline(h = 15000, col = "red", lwd = 2)
dev.off()

pdf(file = paste(PATH_results, "pg-boxplot.pdf", sep=""), width = 4, height = 4)
boxplot(sorted_counts)
dev.off()

#-------------------------------------------------------------------------------

library("org.Mm.eg.db")
library("org.Hs.eg.db")
library(biomaRt)

convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",  host = "https://dec2021.archive.ensembl.org/")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",  host = "https://dec2021.archive.ensembl.org/")
  
  genesV2 <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol",
                    values = x, mart = human,
                    attributesL = c("mgi_symbol"), martL = mouse, uniqueRows = TRUE)
  
  humanx <- unique(genesV2)
  return(humanx)
}

#-------------------------------------------------------------------------------

# Heatmap of neuronal/myelin genes on interest
neurons <- read.csv("./data/neuronal-genes.csv", header = FALSE) #from hDRG prot paper

genelist <- neurons$V1
genelist <- convertHumanGeneList(genelist)

neurons <- trimws(as.character(genelist$MGI.symbol))

data <- as.data.frame(df) %>%
  mutate(GeneID = gsub("\\n", "", rownames(df))) %>%
  distinct(GeneID, .keep_all = TRUE)

data <- data[data$GeneID %in% neurons, ]

head(data)
rownames(data) <- data$GeneID
data$GeneID <- NULL

data <- as.matrix(data)

scaled_expression <- t(scale(t(data), center = TRUE))

match_index <- match(colnames(scaled_expression), colData$sampleID)
colData_reordered <- colData[match_index, ]

tissue_list <- as.factor(colData$Concentration)

# Create a color mapping for colData
tissue_colors <- c("100" = "#3b92df",
                   "25" = "#bed1e1",
                   "50" = "#a2acd9",
                   "Veh" = "#f7e1ae"
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
                                   cluster_rows = TRUE,
                                   cluster_columns = TRUE,
                                   show_row_names = TRUE,
                                   show_column_names = FALSE, #set to TRUE to double check colour legend
                                   row_title = "Proteins",
                                   row_dend_side = "left",
                                   top_annotation = HeatmapAnnotation(tissue = tissue_list, col = list(tissue = col_fun))
)

draw(ht_list, heatmap_legend_side = "right")

#PATH_results = "./output/"

pdf(file = paste0(PATH_results, "/neuronal-heatmap.pdf"), height = 5, width = 6)
draw(ht_list, heatmap_legend_side = "right")
dev.off()

#-------------------------------------------------------------------------------

head(df)

mat <- as.matrix(df[, which(colnames(df) %in% colData$sampleID)])
mat <- mat[complete.cases(mat), ]

rv     <- matrixStats::rowVars(mat) # calculate variance per row (ie. per gene)
select <- order(rv, decreasing=TRUE)[seq_len(min(2000, length(rv)))]

pca <- prcomp(t(mat[select,]), center = TRUE, scale. = TRUE)
pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)

index   <- match(colnames(mat), colData$sampleID)
colData <- colData[index, ]

Concentration   <- as.factor(colData$Concentration)

g <- ggbiplot(pca, choices = c(1,2), 
               groups = interaction(Concentration),
               ellipse = TRUE,
               ellipse.prob = 0.95,
               labels = NULL,
               point.size = 4,
               labels.size = 4, alpha = 1, var.axes = FALSE,
               circle  = TRUE, circle.prob = 0.5,
               varname.size = 3,
               varname.adjust = 1.5,
               varname.abbrev = FALSE)

g <- g + theme(legend.position = 'right') + theme_bw() + coord_fixed(ratio=0.4)

pdf(file = paste(PATH_results, "pca.pdf", sep=""), width = 4, height = 4)
print(g)
dev.off()

#-------------------------------------------------------------------------------

# select distinct genes for testing
data <- as.data.frame(df) %>%
  mutate(GeneID = sub(";.*", "", rownames(df))) %>%  # Remove everything after the semicolon
  distinct(GeneID, .keep_all = TRUE)

rownames(data) <- data$GeneID
data$GeneID <- NULL

mat <- as.matrix(data)

# Verify dimensions of the filtered matrix
dim(mat)

# colData MUST match order of mat
head(colnames(mat))
head(colData$sampleID)

Concentration <- as.factor(colData$Concentration)

design <- model.matrix(~ 0 + Concentration, data = colData) #select design
head(design)

fit <- lmFit(mat, design)

colnames(fit) #check possible comparisons

#----------------

# Oxaliplatin
contrast_matrix <- makeContrasts( Concentration100 -  ConcentrationVeh, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust="BH", number=Inf)

head(results)

sig_de <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 0.5, , drop = FALSE] #extract significant proteins if there
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
  ggtitle("Oxaliplatin")

print(volc)

pdf(file = paste0(PATH_results, "volcano-ox.pdf"), height = 4, width = 4)
print(volc)
dev.off()

pdf(file = paste0(PATH_results, "volcano-ox_big.pdf"), height = 6, width = 6)
print(volc)
dev.off()

table(sig_de$logFC > 1)
write.csv(results, paste(PATH_results, "DEP-analysis-limma_ox.csv"))
write.csv(sig_de,  paste(PATH_results, "DEP-analysis-limma_ox_sig.csv"))

#-------------------------------------------------------------------------------

yang <- read.csv("./data/published/Yang2022_OxaVsVeh_shortterm_stable2.csv")

head(yang)
yang <- yang[, colnames(yang) %in% c("Gene_Symbol", "log2FC", "log2.protein.intensity", "Adjusted_P_value"), ]

head(results)

# # select only DEPs from each
inhouse <- results[results$adj.P.Val < 0.05, ]
yang    <- yang[yang$Adjusted_P_value < 0.05, ]
 
# examine all genes tested
inhouse <- results
yang    <- yang

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
    x = "logFC",
    y = "logFC (Yang2022)"
  ) +
  theme_bw()

print(g)

pdf(file = paste0(PATH_results, "yang2022_correlation_p0.05_100.pdf"), height = 4, width = 4)
print(g)
dev.off()

write.csv(merged_data, "./output/explants/DEPs_across_datasets_withcluster.csv")

clusters <- read.csv("./data/published/Yang2022_Clusters.csv")

merged_data <- merge(merged_data, clusters, by.x = "GeneID", by.y = "Gene.Symbol", all = FALSE)

correlations <- merged_data %>%
  group_by(Condition) %>%
  summarize(correlation = cor(explant.logFC, yang.logFC, use = "complete.obs"), .groups = "drop")

g <- ggplot(merged_data, aes(x = explant.logFC, y = yang.logFC)) 
g <- g + geom_point(color = "#9357f1") + facet_wrap(. ~ Condition)
g <- g + geom_smooth(method = "lm", color = "grey", se = TRUE) 
g <- g + labs(title = paste("DEPs, vs Yang et al 2022.
Corr=",round(correlation, 2)),
              x = "logFC",
              y = "logFC (Yang2022)"
) + theme_bw()
g <- g + geom_text(
  data = correlations,
  aes(label = paste("Corr =", round(correlation, 2))),
  x = Inf, y = -Inf, hjust = 1.1, vjust = -1.1, inherit.aes = FALSE
)

print(g)

pdf(file = paste0(PATH_results, "yang2022_correlation_byCluster.pdf"), height = 4, width = 4)
print(g)
dev.off()

#-------------------------------------------------------------------------------

yang <- read.csv("./data/published/Yang2022_OxaVsVeh_shortterm_stable2.csv")

head(yang)
yang <- yang[, colnames(yang) %in% c("Gene_Symbol", "log2FC", "log2.protein.intensity", "Adjusted_P_value"), ]

head(results)

# examine all genes tested
inhouse <- results
yang    <- yang

merged_data <- merge(inhouse, yang, by.x = "row.names", by.y = "Gene_Symbol", all = FALSE)
head(merged_data)

#intuitive naming for saved csv
colnames(merged_data) <- c("GeneID", "explant.logFC", "explant.AveExpr", "explant.t", 
                           "explant.P.Value", "explant.adj.P.Val", "explant.B",
                           "yang.logFC", "yang.AveExpr", "yang.adj.P.Val")

correlation <- cor(merged_data$explant.AveExpr, merged_data$yang.AveExpr, use = "complete.obs")

g <- ggplot(merged_data, aes(x = explant.AveExpr, y = yang.AveExpr)) 
g <- g + geom_point(color = "#9357f1") 
g <- g + geom_smooth(method = "lm", color = "grey", se = TRUE) 
g <- g + labs(title = paste("DEPs, vs Yang et al 2022.
Corr=",round(correlation, 2)),
              x = "Expression",
              y = "Expression (Yang2022)"
) +
  theme_bw()

print(g)

pdf(file = paste0(PATH_results, "yang2022_correlation_expression.pdf"), height = 4, width = 4)
print(g)
dev.off()

write.csv(merged_data, "./output/explants/Expression_across_datasets.csv")

#-------------------------------------------------------------------------------


