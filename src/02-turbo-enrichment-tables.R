library(dplyr)
library(tidyr)
library(ggplot2)
library("org.Mm.eg.db")
library("org.Hs.eg.db")
library(biomaRt)

PATH_results = "./output/bulk.comparison/"

enrichments   <- read.csv("./output/enrichments.csv", row.names = 1)
bulk_df       <- read.csv("./data/wcl-matrix.csv", row.names = 1)
bulk_colData  <- read.csv("./data/wcl-colData.csv", row.names = 1)
turbo_df      <- read.csv("./data/matrix-for-limma.csv", header = TRUE)
turbo_colData <- read.csv("./data/colData-for-limma.csv", header = TRUE)
  
  convertMouseGeneList <- function(x){
    require("biomaRt")
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",  host = "https://dec2021.archive.ensembl.org/")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",  host = "https://dec2021.archive.ensembl.org/")
    
    genesV2 <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol",
                      values = x, mart = mouse,
                      attributesL = c("hgnc_symbol"), martL = human, uniqueRows = TRUE)
    
    humanx <- unique(genesV2)
    return(humanx)
  }  
  
turbo_df <- turbo_df %>% distinct(genes, .keep_all = TRUE)

drg_df <- turbo_df[, colnames(turbo_df) %in% turbo_colData$sampleID[turbo_colData$Tissue == "DRG"] & 
                     colnames(turbo_df) %in% turbo_colData$sampleID[turbo_colData$Turbo == "T"] | 
                     colnames(turbo_df) == "genes"]


lsc_df <- turbo_df[, colnames(turbo_df) %in% turbo_colData$sampleID[turbo_colData$Tissue == "LSC"] & 
                     colnames(turbo_df) %in% turbo_colData$sampleID[turbo_colData$Turbo == "T"] | 
                     colnames(turbo_df) == "genes"]

paw_df <- turbo_df[, colnames(turbo_df) %in% turbo_colData$sampleID[turbo_colData$Tissue == "paw"] & 
                     colnames(turbo_df) %in% turbo_colData$sampleID[turbo_colData$Turbo == "T"] | 
                     colnames(turbo_df) == "genes"]

scn_df <- turbo_df[, colnames(turbo_df) %in% turbo_colData$sampleID[turbo_colData$Tissue == "SCN"] & 
                     colnames(turbo_df) %in% turbo_colData$sampleID[turbo_colData$Turbo == "T"] | 
                     colnames(turbo_df) == "genes"]


drg_df <- drg_df[drg_df$genes %in% enrichments$Gene[enrichments$Tissue == "DRG"],]
lsc_df <- lsc_df[lsc_df$genes %in% enrichments$Gene[enrichments$Tissue == "LSC"],]
paw_df <- paw_df[paw_df$genes %in% enrichments$Gene[enrichments$Tissue == "paw"],]
scn_df <- scn_df[scn_df$genes %in% enrichments$Gene[enrichments$Tissue == "SCN"],]

library(xlsx)
write.xlsx(drg_df, file="data/enrichment-matrices.xlsx", sheetName="drg", row.names=FALSE)
write.xlsx(lsc_df, file="data/enrichment-matrices.xlsx", sheetName="lsc", append=TRUE, row.names=FALSE)
write.xlsx(paw_df, file="data/enrichment-matrices.xlsx", sheetName="paw", append=TRUE, row.names=FALSE)
write.xlsx(scn_df, file="data/enrichment-matrices.xlsx", sheetName="scn", append=TRUE, row.names=FALSE)

#-------------------------------------------------------------------------------

# combined background across bulk and turbo dfs
background <- read_excel("./data/JRS_curated/20240826_corBackground_for_GO_WCL-TurboALL.xlsx", col_names = FALSE)
names(background)[names(background) == "...1"] <- "genes"

background <- as.data.frame(background)

full_df <- merge(background, drg_df, by = "genes", all = TRUE)
full_df <- merge(full_df, lsc_df, by = "genes", all = TRUE)
full_df <- merge(full_df, paw_df, by = "genes", all = TRUE)
full_df <- merge(full_df, scn_df, by = "genes", all = TRUE)

write.csv(full_df, "./data/turbo_enriched_matrix.csv")

bulk_df <- bulk_df %>% distinct(GeneSymbol, .keep_all = TRUE)
bulk_df$genes <- bulk_df$GeneSymbol
bulk_df$GeneSymbol <- NULL

full_df <- merge(full_df, bulk_df, by = "genes", all.x = TRUE)

symbol.list <- convertMouseGeneList(full_df$genes)
symbol.list$genes <- symbol.list$MGI.symbol

full_df <- merge(full_df, symbol.list, by = "genes", all.x = TRUE)

write.csv(full_df, "./data/full_enriched_matrix.csv")

#-------------------------------------------------------------------------------

# devtools::install_github("lydiaMyr/ImmuCellAI@main")
# BiocManager::install("GSVA")

library(ImmuCellAI)

full_df <- full_df %>% distinct(HGNC.symbol, .keep_all = TRUE)
full_df <- full_df[!is.na(full_df$HGNC.symbol),]

row.names(full_df) <- full_df$HGNC.symbol
full_df$HGNC.symbol <- NULL
full_df$MGI.symbol  <- NULL
full_df$Symbol      <- NULL

full_df[is.na(full_df)] <- 0

immune_enrichment <- ImmuCellAI(sample = full_df, 
           data_type = "rnaseq", 
           group_tag = 1, 
           response_tag = 0, 
           customer = 0,
           min.sz > 1)

abundance_df <- as.data.frame(immune_enrichment$Sample_abundance)

colData <- data.frame(sampleID = bulk_colData$sampleID,
                      Tissue   = bulk_colData$Tissue)

turbo_colData <- turbo_colData[turbo_colData$Turbo == "T", ]
turbo_colData$Tissue <- paste0(turbo_colData$Tissue, "_turbo")

colData2 <- data.frame(sampleID = turbo_colData$sampleID,
                       Tissue   = turbo_colData$Tissue)

colData <- rbind(colData, colData2)

head(colData)

head(abundance_df)
abundance_df$sampleID <- rownames(abundance_df)

df_wide <- merge(colData, abundance_df, by = "sampleID")

g <- ggplot(df_wide, aes(x = Tissue, y = InfiltrationScore)) 
g <- g + geom_boxplot(fill = "skyblue", alpha = 0.7) 
g <- g + geom_jitter(width = 0.2, size = 2, alpha = 0.8) 
g <- g + theme_bw() +
  labs(
    title = " ",
    x = "",
    y = "Immune Cell Infiltration Score"
  ) 
g <- g + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))

print(g)

write.csv(abundance_df, "./output/immunesignature/ImmuCellAI.csv")

PATH_results = "./output/immunesignature/"

pdf(file = paste0(PATH_results, "/ImmuCellAI_boxplot.pdf"), height = 4, width = 4)
print(g)
dev.off()
