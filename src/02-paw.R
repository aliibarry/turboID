library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)


dir.create("./output/paw")
PATH_results = "./output/paw/"

enrichments  <- read.csv("./output/enrichments.csv", row.names = 1)
df      <- read.csv("./data/matrix-for-limma.csv", header = TRUE)
colData <- read.csv("./data/colData-for-limma.csv", header = TRUE)

df <- df %>%
  distinct(genes, .keep_all = TRUE)
rownames(df) <- df$genes

df.paw <- df[, colnames(df) %in% colData$sampleID[colData$Tissue == "LSC"]]
df.paw <- df.paw[, colnames(df.paw) %in% colData$sampleID[colData$Turbo == "T"]]

df.paw <- df.paw[rownames(df.paw) %in% enrichments$Gene[enrichments$Tissue=="LSC"], ]

head(df.paw)

medians_df <- data.frame(Gene = rownames(df.paw), 
                         Median = apply(df.paw, 1, median, na.rm = TRUE))

# write.csv(medians_df, "./output/LSC-proteome-medians.csv")

head(medians_df)


library("org.Mm.eg.db")
library(biomaRt)

genes_of_interest <- read.csv("../nageottes/targetted-protein-list.csv", row.names = 1)
genes_of_interest <- genes_of_interest[genes_of_interest$reason != "cell.type.markers"]

convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",  host = "https://dec2021.archive.ensembl.org/")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",  host = "https://dec2021.archive.ensembl.org/")
  
  #mgi_symbol used instead of ensembl_gene_id for attributes + filters if needed
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x, 
                   mart  = human, attributesL = c("ensembl_gene_id", "mgi_symbol"), 
                   martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2)
  
  return(humanx)
}

hg_genelist <- genes_of_interest$hgnc_symbol
ms_genelist <- convertMouseGeneList(hg_genelist)

#GPCRs/associated 
#ion channels/associated
gois <- ms_genelist$MGI.symbol[ms_genelist$HGNC.symbol %in% genes_of_interest$hgnc_symbol[genes_of_interest$reason == "ion channels/associated"]]
gois2 <- c("Gabra1", "Gabra2", "Gabra5", "Gabrb1", "Gabrb2", "Gabrb3", "Gabrg2", "Gabbr2", "Gabbr1")

gois <- append(gois,gois2)

medians_df$Label <- ifelse(medians_df$Gene %in% gois, medians_df$Gene, NA)
medians_df$LabelColor <- ifelse(!is.na(medians_df$Label), "GPCR", "Other")

# Plot with boxed colored labels
# g1 <- ggplot(medians_df, aes(x = reorder(Gene, Median), y = Median)) +
#   geom_point(color = "grey40") +
#   geom_label_repel(
#     aes(label = Label, fill = LabelColor),
#     color = "black",  # text color
#     box.padding = 0.3,
#     max.overlaps = Inf,
#     show.legend = FALSE
#   ) +
#   scale_fill_manual(values = c("GPCR" = "#98c1d6", "Other" = NA)) +
#   labs(
#     title = "Central Terminals",
#     x = "Rank",
#     y = "Median Expression"
#   ) +
#   theme_minimal() +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

g2 <- ggplot(medians_df, aes(x = reorder(Gene, Median), y = Median)) +
  geom_point(color = "grey40") +
  geom_label_repel(
    aes(label = Label, fill = LabelColor),
    color = "black",  # text color
    box.padding = 0.3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("GPCR" = "#b1e286", "Other" = NA)) +
  labs(
    title = "Central Terminals",
    x = "Rank",
    y = "Median Expression"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

print(g1)
print(g2)

library(patchwork)

PATH_results = "../nageottes/output/"

pdf(file = paste0(PATH_results, "/TurboID-centralterminals.pdf"), height = 8, width = 8)
print(g1/g2)
dev.off()

