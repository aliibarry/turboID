library(dplyr)
library(ggplot2)

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


# Extract ion channels from published list of receptors
goi <- read.table("../hdrg_proteomics/ionchannels.csv", fill = TRUE, header = TRUE, sep=",")

genelist <- goi$HGNC.symbol #extract symbols
genelist <- convertHumanGeneList(genelist)

goi <- merge(genelist, goi, 
             by.x = "HGNC.symbol", 
             by.y = "HGNC.symbol", 
             all.x = TRUE)
goi <- goi %>% distinct() # remove duplicates
head(goi)

ion.channels <- goi %>% 
  filter(str_detect(MGI.symbol, "^Trp") | 
           str_detect(MGI.symbol, "^Kcn") | 
           str_detect(MGI.symbol, "^Scn") |
           str_detect(MGI.symbol, "^Gab") | 
           str_detect(MGI.symbol, "^Asic") | 
           str_detect(MGI.symbol, "^Cng") | 
           str_detect(MGI.symbol, "^Ryr") | 
           str_detect(MGI.symbol, "^Cac") | 
           str_detect(MGI.symbol, "^Hvc"))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Cross-tissue dataset
df           <- read.csv("./data/matrix-for-limma.csv", row.names = 1, check.names = FALSE)
colData      <- read.csv("./data/colData-for-limma.csv", row.names = 1)
enrichments  <- read.csv("./output/enrichments_75filt.csv", row.names = 1)

colData <- colData[colData$Turbo == "T", ] # optional, change plot output name to match

# filter and reorder
mat <- df[rownames(df) %in% enrichments$Gene, ]
mat$genes <- NULL

mat     <- mat[,colnames(mat) %in% colData$sampleID]
colData <- colData[colData$sampleID %in% colnames(mat), ]
index   <- match(colnames(mat), colData$sampleID)
colData <- colData[index, ]

head(colnames(mat))
head(colData$sampleID)

# trim gene names for whitespace to ensure matches
rownames(mat) <- trimws(as.character(rownames(mat)))

data <- mat[rownames(mat) %in% ion.channels$MGI.symbol, ]

# prep data for heatmap
scaled_expression <- t(scale(t(data), center = TRUE))

match_index <- match(colnames(scaled_expression), colData$sampleID)
colData_reordered <- colData[match_index, ]

tissue_list <- as.factor(paste(colData_reordered$Tissue,"-",colData_reordered$Turbo))

# # Create a color mapping for metadata
# tissue_colors <- c("DRG - T" = "#f55c3a",
#                    "DRG - TC" = "#f1ccc4",
#                    "LSC - T" = "#edb127",
#                    "LSC - TC" = "#f7e1ae",
#                    "SCN - T" = "#29c99d",
#                    "SCN - TC" = "#abcfc5",
#                    "paw - T" = "#95459b",
#                    "paw - TC" = "#e7d9e8"
# )

# with Turbo controls removed
tissue_colors <- c("DRG - T" = "#f55c3a",
                   "LSC - T" = "#edb127",
                   "SCN - T" = "#29c99d",
                   "paw - T" = "#95459b"
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

PATH_results = "./output/"

pdf(file = paste0(PATH_results, "ionchannel-heatmap_Turbo.pdf"), height = 7, width = 8)
draw(ht_list, heatmap_legend_side = "right")
dev.off()

#-------------------------------------------------------------------------------

# Explant data plots
df           <- read.csv("./data/explants/combined_genes.csv", row.names = 1, check.names = FALSE)
colData      <- read.csv("./data/explants/colData.csv", row.names = 1)
enrichments  <- read.csv("./output/explants/tissue-enr/enrichments.csv", row.names = 1)

colData <- colData[colData$Turbo == "T", ] # optional, change plot output name to match

mat <- df[rownames(df) %in% enrichments$Gene, ]
mat$genes <- NULL

mat     <- mat[,colnames(mat) %in% colData$sampleID]
colData <- colData[colData$sampleID %in% colnames(mat), ]
index   <- match(colnames(mat), colData$sampleID)
colData <- colData[index, ]

head(colnames(mat))
head(colData$sampleID)

rownames(mat) <- trimws(as.character(rownames(mat)))

data <- mat[rownames(mat) %in% ion.channels$MGI.symbol, ]

scaled_expression <- t(scale(t(data), center = TRUE))

match_index <- match(colnames(scaled_expression), colData$sampleID)
colData_reordered <- colData[match_index, ]

tissue_list <- as.factor(paste(colData_reordered$Condition,"-",colData_reordered$Turbo))

# # Create a color mapping for metadata
# tissue_colors <- c("Ox - T" = "#3b92df",
#                    "Ox - TC" = "#bed1e1",
#                    "V - T" = "#edb127",
#                    "V - TC" = "#f7e1ae"
# )

tissue_colors <- c("Ox - T" = "#3b92df",
                   "V - T" = "#edb127"
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

# Create Heatmap
ht_list <- ComplexHeatmap::Heatmap(scaled_expression,
                                   #name = "Expression",
                                   col= col_fun2,
                                   clustering_distance_columns = "manhattan",
                                   cluster_rows = TRUE,
                                   cluster_columns = TRUE,
                                   show_row_names = TRUE,
                                   show_column_names = TRUE, #set to TRUE to double check colour legend
                                   row_title = "Proteins",
                                   row_dend_side = "left",
                                   top_annotation = HeatmapAnnotation(tissue = tissue_list, col = list(tissue = col_fun))
)

draw(ht_list, heatmap_legend_side = "right")

PATH_results = "./output/explants/"

pdf(file = paste0(PATH_results, "ionchannel-heatmap_Turbo.pdf"), height = 5, width = 6)
draw(ht_list, heatmap_legend_side = "right")
dev.off()
