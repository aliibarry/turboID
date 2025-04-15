library(dplyr)
library(tidyr)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)

file_path1 <- "data/published/Paw.SCN.Supplementary Table1.xlsx" #STable 1 doi: 10.7554/eLife.81431, paw
file_path2 <- "data/published/Paw.SCN.Supplementary Table2.xlsx" #STable 2 doi: 10.7554/eLife.81431, scn
sheet1.1 <- read_excel(file_path1, sheet = 1) # sheets 2-8 for DEPs
sheet1.2 <- read_excel(file_path2, sheet = 1) # sheets 2-8 for DEPs

df      <- read.csv("./data/matrix-for-limma.csv", header = TRUE)
colData <- read.csv("./data/colData-for-limma.csv", header = TRUE)

PATH_results = "./output/published.comparison/"

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Remove the first 2 rows
paw_df <- sheet1.1[-c(1, 2), ]

# Display the cleaned data
head(paw_df)

colnames(paw_df) <- as.character(paw_df[1, ])
paw_df <- paw_df[-1, ]  # Remove the row that is now the header
paw_df <- as.data.frame(paw_df[, -ncol(paw_df)])

head(paw_df)

paw_df <- paw_df %>%
  separate(`Gene Name`, into = paste0("genes"), sep = ";", remove = TRUE)

paw_df$genes <- trimws(as.character(paw_df$genes))

paw_df <- paw_df %>%
  mutate(genes = gsub("\\n", "", genes)) %>%
  distinct(genes, .keep_all = TRUE)

head(paw_df)

paw_genes <- paw_df$genes

#-------------------------------------------------------------------------------

df <- df %>%
  mutate(genes = gsub("\\n", "", genes)) %>%
  distinct(genes, .keep_all = TRUE)

rownames(df) <- df$genes

paw_samples <- colData[colData$Tissue == 'paw', ]
paw_data <- df[, paw_samples$sampleID, drop = FALSE]
detected_genes_paw <- rownames(paw_data)[rowSums(!is.na(paw_data)) > 0]

common_genes <- intersect(detected_genes_paw, paw_genes)
unique_genes_turbo <- setdiff(detected_genes_paw, common_genes)


enrichment_paw <- enrichments$Gene[enrichments$Tissue == "paw"]
enriched_common <- intersect(paw_genes, enrichment_paw)
enriched_turbo  <- setdiff(enrichment_paw, common_genes)

write.csv(enriched_turbo, "./output/turbo-specpfic-vs-published.csv")


# plot_data <- data.frame(
#   Category = c("WT.published", "Overlap", "TurboID"),
#   Count = c(
#     length(paw_genes) - length(common_genes),  # Count of common genes
#     length(common_genes),  # Overlap (common genes)
#     length(detected_genes_paw) - length(common_genes)  # Specific to detected genes in paw
#   ),
#   Tissue = "paw"
# )
# 
# # Reshape the data to fit ggplot2 requirements
# plot_data <- plot_data %>%
#   mutate(Category = factor(Category, levels = c("TurboID", "Overlap", "WT.published")))

plot_data <- data.frame(
  Category = c("Overlap", "TurboID"),
  Count = c(
    length(common_genes),  # Overlap (common genes)
    length(detected_genes_paw) - length(common_genes)  # Specific to detected genes in paw
  ),
  Tissue = "paw"
)

# Reshape the data to fit ggplot2 requirements
plot_data <- plot_data %>%
  mutate(Category = factor(Category, levels = c("TurboID", "Overlap")))

# Create the stacked bar plot
g <- ggplot(plot_data, aes(x = Tissue, y = Count, fill = Category)) 
g <- g + geom_bar(stat = "identity", position = "stack") 
g <- g + geom_text(
    aes(label = Count),
    position = position_stack(vjust = 0.5),
    color = "black"
  ) 
g <- g + labs(title = "", x = "", y = "PG Count") +
  theme_bw() +
  scale_fill_manual(
    values = c("WT.published" = "#f6c8f4ff", "Overlap" = "#c677c2ff", "TurboID" = "#6b3268ff"),
    name = "Category"
  )

# Print the plot
print(g)

pdf(file = paste(PATH_results, "PG_counts_paw.enriched_noWT.pdf", sep=""), width = 4, height = 4)
print(g)
dev.off()


paw_turbo <- data.frame(GeneID = unique_genes_turbo)
# write.csv(paw_turbo, paste0(PATH_results, "paw.turbo.csv"))

enrichments  <- read.csv("./output/enrichments.csv", row.names = 1)

paw_turbo_enriched <- paw_turbo[paw_turbo$GeneID %in% enrichments$Gene, ]
 

#-------------------------------------------------------------------------------


background <- rbind(data.frame(genes = unique_genes_turbo), data.frame(genes = paw_genes))

ego <- enrichGO(gene          = unique_genes_turbo,
                universe      = background$genes,
                OrgDb         = org.Mm.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
#head(ego)
goplot(ego)
enrichMap(ego)

pdf(file = paste(PATH_results, "BP-turbo.paw-network.pdf", sep=""), width = 8, height = 8)
goplot(ego)
dev.off()

pdf(file = paste(PATH_results, "BP-turbo.paw-barplot.pdf", sep=""), width = 6, height = 3)
mutate(ego, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")
dev.off()

write.csv(ego, paste(PATH_results, "BP-turbo.paw.csv"))


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Extract SCN data for comparison

# Remove the first 2 rows
scn_df <- sheet1.2[-c(1, 2), ]

# Display the cleaned data
head(scn_df)

colnames(scn_df) <- as.character(scn_df[1, ])
scn_df <- scn_df[-1, ]  # Remove the row that is now the header
#scn_df <- as.data.frame(scn_df[, -ncol(scn_df)])

head(scn_df)

scn_df <- scn_df %>%
  separate(`Gene Name`, into = paste0("genes"), sep = ";", remove = TRUE)

scn_df$genes <- trimws(as.character(scn_df$genes))

scn_df <- scn_df %>%
  mutate(genes = gsub("\\n", "", genes)) %>%
  distinct(genes, .keep_all = TRUE)

head(scn_df)

scn_genes <- scn_df$genes

#-----------

df <- df %>%
  mutate(genes = gsub("\\n", "", genes)) %>%
  distinct(genes, .keep_all = TRUE)

rownames(df) <- df$genes

scn_samples <- colData[colData$Tissue == 'SCN', ]
scn_data <- df[, scn_samples$sampleID, drop = FALSE]
detected_genes_scn <- rownames(scn_data)[rowSums(!is.na(scn_data)) > 0]

common_genes <- intersect(detected_genes_scn, scn_genes)
unique_genes_turbo <- setdiff(detected_genes_scn, common_genes)

# plot_data <- data.frame(
#   Category = c("WT.published", "Overlap", "TurboID"),
#   Count = c(
#     length(scn_genes) - length(common_genes),  # Count of common genes
#     length(common_genes),  # Overlap (common genes)
#     length(detected_genes_scn) - length(common_genes)  # Specific to detected genes in paw
#   ),
#   Tissue = "SCN"
# )
# 
# # Reshape the data to fit ggplot2 requirements
# plot_data <- plot_data %>%
#   mutate(Category = factor(Category, levels = c("TurboID", "Overlap", "WT.published")))


plot_data <- data.frame(
  Category = c("Overlap", "TurboID"),
  Count = c(
    length(common_genes),  # Overlap (common genes)
    length(detected_genes_scn) - length(common_genes)  # Specific to detected genes in paw
  ),
  Tissue = "SCN"
)

# Reshape the data to fit ggplot2 requirements
plot_data <- plot_data %>%
  mutate(Category = factor(Category, levels = c("TurboID", "Overlap")))

# Create the stacked bar plot
g <- ggplot(plot_data, aes(x = Tissue, y = Count, fill = Category)) 
g <- g + geom_bar(stat = "identity", position = "stack") 
g <- g + geom_text(
  aes(label = Count),
  position = position_stack(vjust = 0.5),
  color = "black"
) 
g <- g + labs(title = "", x = "", y = "PG Count") +
  theme_bw() +
  scale_fill_manual(
    values = c("WT.published" = "#d1dbfdff", "Overlap" = "#899acfff", "TurboID" = "#455ba3ff"),
    name = "Category"
  )

# Print the plot
print(g)

pdf(file = paste(PATH_results, "PG_counts_scn.enriched_noWT.pdf", sep=""), width = 4, height = 4)
print(g)
dev.off()

#-----------

background <- rbind(data.frame(genes = unique_genes_turbo), data.frame(genes = scn_genes))

ego <- enrichGO(gene          = unique_genes_turbo,
                universe      = background$genes,
                OrgDb         = org.Mm.eg.db,
                keyType       = "SYMBOL",
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
#head(ego)
goplot(ego)
enrichMap(ego)

pdf(file = paste(PATH_results, "CC-turbo.scn-network.pdf", sep=""), width = 8, height = 8)
goplot(ego)
dev.off()

pdf(file = paste(PATH_results, "CC-turbo.scn-barplot.pdf", sep=""), width = 6, height = 3)
mutate(ego, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")
dev.off()

write.csv(ego, paste(PATH_results, "CC-turbo.scn.csv"))

scn_turbo <- data.frame(GeneID = unique_genes_turbo)
write.csv(scn_turbo, paste0(PATH_results, "SCN.turbo.csv"))









