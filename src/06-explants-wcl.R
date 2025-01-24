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

write.csv(df, "./data/explants/wcl/matrix.csv")
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
                         fontsize_col = 8
)

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

# data <- df %>%
#   separate(GeneID, into = paste0("GeneID"), sep = ";", remove = TRUE)

df$GeneID <- trimws(as.character(df$GeneID))

data <- df %>%
  mutate(GeneID = gsub("\\n", "", GeneID)) %>%
  distinct(GeneID, .keep_all = TRUE)

# data$genes <- NULL
# expression_data <- data[complete.cases(data),]

data <- data[data$GeneID %in% neurons, ]

head(data)
rownames(data) <- data$GeneID
data$GeneID <- NULL

data <- as.matrix(data)

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


