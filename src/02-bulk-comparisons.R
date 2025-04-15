library(dplyr)
library(tidyr)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)


PATH_results = "./output/bulk.comparison/"

enrichments  <- read.csv("./output/enrichments.csv", row.names = 1)
bulk_df      <- read.csv("./data/wcl-matrix.csv", row.names = 1, check.names = FALSE)
bulk_colData <- read.csv("./data/wcl-colData.csv", row.names = 1)

#-------------------------------------------------------------------------------

# get specific geneIDs
bulk_df <- bulk_df %>%
  separate(`GeneID`, into = paste0("GeneID"), sep = ";", remove = TRUE)

bulk_df <- bulk_df %>%
  mutate(GeneID = gsub("\\n", "", GeneID)) %>%
  distinct(GeneID, .keep_all = TRUE)

rownames(bulk_df) <- bulk_df$GeneID
bulk_df$GeneID <- NULL

#-------------------------------------------------------------------------------

# check for unique genes in paw turbo sample (vs bulk proteomics)
paw_samples <- bulk_colData[bulk_colData$Tissue == 'paw', ]
paw_data    <- bulk_df[, paw_samples$sampleID, drop = FALSE]
paw_data    <- paw_data[rowSums(!is.na(paw_data)) > 0, ]
detected_genes <- rownames(paw_data)[rowSums(!is.na(paw_data)) > 0]

paw_genes <- enrichments$Gene[enrichments$Tissue == "paw"]

common_genes <- intersect(detected_genes, paw_genes)
unique_turbo <- setdiff(paw_genes, common_genes)

plot_data <- data.frame(
  Category = c("Overlap", "TurboID"),
  Count = c(
     length(common_genes),
     length(paw_genes) - length(common_genes)
  ),
  Tissue = "paw"
)

plot_data <- plot_data %>%
  mutate(Category = factor(Category, levels = c("TurboID", "Overlap")))

# Create the stacked bar plot
g <- ggplot(plot_data, aes(x = Tissue, y = Count, fill = Category)) 
g <- g + geom_bar(stat = "identity", position = "stack") 
g <- g + geom_text(
    aes(label = Count),
    position = position_stack(vjust = 0.5),
    color = "black") 
g <- g + labs(title = "", x = "", y = "PG Count") +
  theme_bw() +
  scale_fill_manual(
    values = c("WT.bulk" = "#f6c8f4ff", "Overlap" = "#c677c2ff", "TurboID" = "#6b3268ff"),
    name = "Category")

print(g)

pdf(file = paste(PATH_results, "PG_counts_paw.enriched.pdf", sep=""), width = 4, height = 4)
print(g)
dev.off()

#-------------------------------------------------------------------------------

# Check TurboID-specific proteins across all tissues

all_plot_data <- list()
tissues <- unique(bulk_colData$Tissue)

for (tissue in tissues) {
  tissue_samples <- bulk_colData[bulk_colData$Tissue == tissue, ]
  tissue_data    <- bulk_df[, tissue_samples$sampleID, drop = FALSE]
  detected_genes <- rownames(tissue_data)[rowSums(!is.na(tissue_data)) > 0]
  tissue_genes   <- enrichments$Gene[enrichments$Tissue == tissue]
  common_genes   <- intersect(detected_genes, tissue_genes)
  unique_turbo   <- setdiff(detected_genes, common_genes)
  
  # Create the plot data for the current tissue
  plot_data <- data.frame(
    Category = c("Overlap", "TurboID"),
    Count = c(
      length(common_genes),
      length(tissue_genes) - length(common_genes)
    ),
    Tissue = tissue
  )
  
  all_plot_data[[tissue]] <- plot_data
}

final_plot_data <- bind_rows(all_plot_data)
final_plot_data <- final_plot_data %>%
  mutate(Category = factor(Category, levels = c("TurboID", "Overlap")))

g <- ggplot(final_plot_data, aes(x = Tissue, y = Count, fill = Category)) 
g <- g +geom_bar(stat = "identity", position = "stack") 
g <- g +  geom_text(
    aes(label = Count),
    position = position_stack(vjust = 0.5),
    color = "black") 
g <- g +labs(title = "", x = "", y = "PG Count") 
g <- g + theme_bw() +
  scale_fill_manual(
    values = c("Overlap" = "#dededeff", "TurboID" = "#5db7b9ff"),
    name = "Category")

# Print the plot
print(g)

pdf(file = paste(PATH_results, "PG_counts_vs_Bulk.pdf", sep=""), width = 5, height = 4)
print(g)
dev.off()

#-------------

# Get geneIDs specific to neuronal-enriched proteomes (vs bulk)
gene_info_list <- list()
for (tissue in tissues) {
  #extract data per tissue and get all identified genes
  tissue_samples <- bulk_colData[bulk_colData$Tissue == tissue, ]
  tissue_data    <- bulk_df[, tissue_samples$sampleID, drop = FALSE]
  detected_genes <- rownames(tissue_data)[rowSums(!is.na(tissue_data)) > 0]
  tissue_genes   <- enrichments$Gene[enrichments$Tissue == tissue]
  
  # genes present in both bulk and neuronally-enriched proteomes
  common_genes <- intersect(detected_genes, tissue_genes)
  
  # turbo-specific
  unique_genes <- setdiff(tissue_genes, common_genes)

  common_df <- data.frame(
    Gene = common_genes,
    Tissue = tissue,
    `Common/Unique` = "Common")

  unique_df <- data.frame(
    Gene = unique_genes,
    Tissue = tissue,
    `Common/Unique` = "Unique")

  gene_info_list[[tissue]] <- rbind(common_df, unique_df)
}

bulk.overlap <- bind_rows(gene_info_list)

uni_df <- bulk.overlap[bulk.overlap$Common.Unique == "Unique", ]
table(uni_df$Tissue)

write_csv(bulk.overlap, paste0(PATH_results, "overlap-with-bulk.csv"))

#-------------------------------------------------------------------------------

#Pathway analyses
background  <- bulk.overlap$Gene[bulk.overlap$Tissue == "DRG"]
unique_list <- bulk.overlap$Gene[bulk.overlap$Tissue == "DRG" & bulk.overlap$Common.Unique == "Unique"]

ego <- enrichGO(gene          = unique_list,
                universe      = background,
                OrgDb         = org.Mm.eg.db,
                keyType       = "SYMBOL",
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
#head(ego)
goplot(ego)







