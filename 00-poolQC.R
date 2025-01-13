library(diann)
library(dplyr)
library(ggplot2)

df <- diann_load("./data/DIANN_output/SCN_pool/SCN_pool_report.tsv")

precursors <- diann_matrix(df, q = 0.01)


peptides.maxlfq <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,], 
                                sample.header = "Run",
                                group.header="Stripped.Sequence", 
                                id.header = "Precursor.Id", 
                                quantity.header = "Precursor.Normalised")

peptides.maxlfq <- as.data.frame(peptides.maxlfq)
peptides.maxlfq$ModifiedSequence <- rownames(peptides.maxlfq)
rownames(peptides.maxlfq) <- NULL
peptides.maxlfq <- peptides.maxlfq %>%
  select(ModifiedSequence, everything())
# head(peptides.maxlfq)

gene.groups <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$GG.Q.Value <= 0.01,], 
                            sample.header = "Run",
                            group.header="Genes", 
                            id.header = "Precursor.Id", 
                            quantity.header = "Precursor.Normalised")

write.csv(precursors, "./data/SCN_precursors.csv", na="NA",eol = "\n", row.names = T)
write.csv(peptides.maxlfq, "./data/SCN_peptides.csv", na="NA", eol = "\n", row.names = FALSE)
write.csv(gene.groups, "./data/SCN_pgs.csv", na="NA", eol = "\n", row.names = T)

#-------------------------------------------------------------------------------

df <- diann_load("./data/DIANN_output/LDRG_pool/LDRG_pool_report.tsv")

precursors <- diann_matrix(df, q = 0.01)


peptides.maxlfq <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,], 
                                sample.header = "Run",
                                group.header="Stripped.Sequence", 
                                id.header = "Precursor.Id", 
                                quantity.header = "Precursor.Normalised")

peptides.maxlfq <- as.data.frame(peptides.maxlfq)
peptides.maxlfq$ModifiedSequence <- rownames(peptides.maxlfq)
rownames(peptides.maxlfq) <- NULL
peptides.maxlfq <- peptides.maxlfq %>%
  select(ModifiedSequence, everything())
# head(peptides.maxlfq)

gene.groups <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$GG.Q.Value <= 0.01,], 
                            sample.header = "Run",
                            group.header="Genes", 
                            id.header = "Precursor.Id", 
                            quantity.header = "Precursor.Normalised")

write.csv(precursors, "./data/DRG_precursors.csv", na="NA",eol = "\n", row.names = T)
write.csv(peptides.maxlfq, "./data/DRG_peptides.csv", na="NA", eol = "\n", row.names = FALSE)
write.csv(gene.groups, "./data/DRG_pgs.csv", na="NA", eol = "\n", row.names = T)

#-------------------------------------------------------------------------------

df <- diann_load("./data/DIANN_output/dLSC_pool/dLSC_pool_report.tsv")

precursors <- diann_matrix(df, q = 0.01)


peptides.maxlfq <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,], 
                                sample.header = "Run",
                                group.header="Stripped.Sequence", 
                                id.header = "Precursor.Id", 
                                quantity.header = "Precursor.Normalised")

peptides.maxlfq <- as.data.frame(peptides.maxlfq)
peptides.maxlfq$ModifiedSequence <- rownames(peptides.maxlfq)
rownames(peptides.maxlfq) <- NULL
peptides.maxlfq <- peptides.maxlfq %>%
  select(ModifiedSequence, everything())
# head(peptides.maxlfq)

gene.groups <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$GG.Q.Value <= 0.01,], 
                            sample.header = "Run",
                            group.header="Genes", 
                            id.header = "Precursor.Id", 
                            quantity.header = "Precursor.Normalised")

write.csv(precursors, "./data/LSC_precursors.csv", na="NA",eol = "\n", row.names = T)
write.csv(peptides.maxlfq, "./data/LSC_peptides.csv", na="NA", eol = "\n", row.names = FALSE)
write.csv(gene.groups, "./data/LSC_pgs.csv", na="NA", eol = "\n", row.names = T)

#-------------------------------------------------------------------------------

df <- diann_load("./data/DIANN_output/paw_pool/paw_pool_report.tsv")

precursors <- diann_matrix(df, q = 0.01)

peptides.maxlfq <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,], 
                                sample.header = "Run",
                                group.header="Stripped.Sequence", 
                                id.header = "Precursor.Id", 
                                quantity.header = "Precursor.Normalised")

peptides.maxlfq <- as.data.frame(peptides.maxlfq)
peptides.maxlfq$ModifiedSequence <- rownames(peptides.maxlfq)
rownames(peptides.maxlfq) <- NULL
peptides.maxlfq <- peptides.maxlfq %>%
  select(ModifiedSequence, everything())
# head(peptides.maxlfq)

gene.groups <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$GG.Q.Value <= 0.01,], 
                            sample.header = "Run",
                            group.header="Genes", 
                            id.header = "Precursor.Id", 
                            quantity.header = "Precursor.Normalised")

write.csv(precursors, "./data/paw_precursors.csv", na="NA",eol = "\n", row.names = T)
write.csv(peptides.maxlfq, "./data/paw_peptides.csv", na="NA", eol = "\n", row.names = FALSE)
write.csv(gene.groups, "./data/paw_pgs.csv", na="NA", eol = "\n", row.names = T)

#-------------------------------------------------------------------------------

PATH_results = "./output/QC/"

df1 <- read.csv("./data/paw_peptides.csv")
df2 <- read.csv("./data/LSC_peptides.csv")
df3 <- read.csv("./data/DRG_peptides.csv")
df4 <- read.csv("./data/SCN_peptides.csv")

peptides <- Reduce(function(x, y) merge(x, y, by = "ModifiedSequence", all = TRUE), list(df1, df2, df3, df4))
str(peptides)

peptides$ModifiedSequence <- NULL

pre_counts  <- colSums(!is.na(peptides))
pre_counts  <- as.data.frame(pre_counts)

sorted_indices <- order(-pre_counts$pre_counts)
sorted_counts  <- pre_counts$pre_counts[sorted_indices]
sorted_names   <- rownames(pre_counts)[sorted_indices]

pdf(file = paste(PATH_results, "pool/counts-precursor.pdf", sep=""), width = 5, height = 7)
par(mfrow=c(2,1))
barplot(sorted_counts, 
        names.arg = sorted_names,
        ylab = "Precursor Count", 
        las = 3,
        cex.names = 0.8)
dev.off()

#-------------------------------------------------------------------------------

df1 <- read.csv("./data/paw_pgs.csv")
df2 <- read.csv("./data/LSC_pgs.csv")
df3 <- read.csv("./data/DRG_pgs.csv")
df4 <- read.csv("./data/SCN_pgs.csv")

pgs <- Reduce(function(x, y) merge(x, y, by = "X", all = TRUE), list(df1, df2, df3, df4))
str(pgs)

pgs$X <- NULL

pre_counts  <- colSums(!is.na(pgs))
pre_counts  <- as.data.frame(pre_counts)

sorted_indices <- order(-pre_counts$pre_counts)
sorted_counts  <- pre_counts$pre_counts[sorted_indices]
sorted_names   <- rownames(pre_counts)[sorted_indices]

pdf(file = paste(PATH_results, "pool/counts-pgs.pdf", sep=""), width = 5, height = 7)
par(mfrow=c(2,1))
barplot(sorted_counts, 
        names.arg = sorted_names,
        ylab = "PG Count", 
        las = 3,
        cex.names = 0.8)
dev.off()

#-------------------------------------------------------------------------------

mat <- pgs[complete.cases(pgs), ]
correlation_matrix <- cor(mat)

# Melt the correlation matrix into long format
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


pdf(file = paste(PATH_results, "pool/correlation.pdf", sep=""), width = 6, height = 5.5)
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
