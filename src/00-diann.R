library(dplyr)
library(tidyr)
library(diann)

#-------------------------------------------------------------------------------

df <- diann_load("/mnt/CommonStorage/AGSchmidt/8_Transfer/20250114_TurboID_DRGexpl_Ox/Ox_T/Ox_T_report.tsv")

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

write.csv(precursors,      "./data/explants/Ox_T_precursors.csv", na="NA",eol = "\n", row.names = T)
write.csv(peptides.maxlfq, "./data/explants/Ox_T_peptides.csv", na="NA", eol = "\n", row.names = T)
write.csv(gene.groups,     "./data/explants/Ox_T_genegroups.csv", na="NA", eol = "\n", row.names = T)

#-------------------------------------------------------------------------------

df <- diann_load("/mnt/CommonStorage/AGSchmidt/8_Transfer/20250114_TurboID_DRGexpl_Ox/Ox_TC/Ox_TC_report.tsv")

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

write.csv(precursors,      "./data/explants/Ox_TC_precursors.csv", na="NA",eol = "\n", row.names = T)
write.csv(peptides.maxlfq, "./data/explants/Ox_TC_peptides.csv", na="NA", eol = "\n", row.names = T)
write.csv(gene.groups,     "./data/explants/Ox_TC_genegroups.csv", na="NA", eol = "\n", row.names = T)

#-------------------------------------------------------------------------------

df <- diann_load("/mnt/CommonStorage/AGSchmidt/8_Transfer/20250114_TurboID_DRGexpl_Ox/V_T/V_T_report.tsv")

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

write.csv(precursors,      "./data/explants/V_T_precursors.csv", na="NA",eol = "\n", row.names = T)
write.csv(peptides.maxlfq, "./data/explants/V_T_peptides.csv", na="NA", eol = "\n", row.names = T)
write.csv(gene.groups,     "./data/explants/V_T_genegroups.csv", na="NA", eol = "\n", row.names = T)

#-------------------------------------------------------------------------------

df <- diann_load("/mnt/CommonStorage/AGSchmidt/8_Transfer/20250114_TurboID_DRGexpl_Ox/Pool_T-TC/Ox_pool_report.tsv")

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

write.csv(precursors,      "./data/explants/pool_precursors.csv", na="NA",eol = "\n", row.names = T)
write.csv(peptides.maxlfq, "./data/explants/pool_peptides.csv", na="NA", eol = "\n", row.names = T)
write.csv(gene.groups,     "./data/explants/pool_genegroups.csv", na="NA", eol = "\n", row.names = T)

#-------------------------------------------------------------------------------

PATH_results = "./output/explants/QC/"

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

genegroups_files <- list.files(path = "./data/explants", pattern = "genegroups.csv", full.names = TRUE)

data_list <- list()
data_list <- lapply(genegroups_files, function(file) {
  data <- read.csv(file, row.names = 1, check.names = FALSE)  # Read file and set row names
  return(data)
})

merged_data <- Reduce(function(x, y) {
  merged <- merge(x, y, by = "row.names", all = TRUE)
  merged <- merged[, !duplicated(colnames(merged))]
  rownames(merged) <- merged$Row.names
  merged <- merged[, -1]  # Remove the redundant Row.names column
  
  return(merged)
}, data_list)

head(merged_data)
df <- log2(merged_data)

# Build corresponding metadata
sample_names  <- colnames(merged_data)
sample_info   <- strsplit(sample_names, "_")
colData       <- as.data.frame(do.call(rbind, sample_info))
colnames(colData) <- c("Condition", "Sex", "Turbo", "Slot", "Value", "Sample")
colData$Turbo <- gsub("[0-9]", "", as.character(colData$Turbo))

colData <- as.data.frame(cbind(sampleID = sample_names, colData))
head(colData)

write.csv(df, "./data/explants/combined_genes.csv")
write.csv(colData, "./data/explants/colData.csv")





