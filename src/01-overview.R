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

dir.create("./output/QC")
PATH_results = "./output/QC/"

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# extract all sheets from an excel document
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

overview   <- read_excel_allsheets("./data/JRS_curated/20240527overview_TurboID_March24.xlsx")
background <- read_excel_allsheets("./data/JRS_curated/20240826_corBackground_for_GO_WCL-TurboALL.xlsx")
paw_rerun  <- read_excel_allsheets("./data/JRS_curated/20241219_rerun_paw_DIANNquant.xlsx")

#-------------------------------------------------------------------------------

df <- read.csv("./data/DIANN_R_output/allTissues_T-TC_Com_gg.csv", check.names = FALSE)

column_names <- colnames(df)[-1]
metadata <- do.call(rbind, strsplit(column_names, "_"))

metadata <- as.data.frame(metadata)
colnames(metadata) <- c("Sample", "Tissue", "Turbo", "Sex", "ID") #clarify TC/T and ID from Julia
metadata$sampleID <- column_names

#-------------------------------------------------------------------------------

# Reshape the dataframe to long format
df_long <- df %>%
  pivot_longer(cols = -GeneID, names_to = "sampleID", values_to = "Expression") %>%
  left_join(metadata, by = "sampleID")

summed_expression <- df_long %>%
  group_by(sampleID, Tissue, Turbo) %>%
  summarize(SumCount = sum(Expression, na.rm = TRUE), .groups = "drop") %>%
  mutate(Tissue = factor(Tissue, levels = c("LSC", "DRG", "SCN", "paw")))

mean_intensity <- df_long %>%
  group_by(sampleID, Tissue, Turbo) %>%
  summarize(MeanIntensity = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
  mutate(Tissue = factor(Tissue, levels = c("LSC", "DRG", "SCN", "paw")))

detected_counts <- df_long %>%
  group_by(sampleID, Tissue, Turbo) %>%
  summarize(DetectedCount = sum(!is.na(Expression)), .groups = "drop") %>%
  mutate(Tissue = factor(Tissue, levels = c("LSC", "DRG", "SCN", "paw")))

#-------

g1 <- ggplot(summed_expression, aes(x = Tissue, y = SumCount, fill = Turbo)) 
g1 <- g1 + geom_boxplot() 
g1 <- g1 + theme_classic() +
  labs(title = "", 
       x = "", 
       y = "Summed Expression")

print(g1)

#------

g2 <- ggplot(mean_intensity, aes(x = Tissue, y = MeanIntensity, fill = Turbo)) 
g2 <- g2 + geom_boxplot() 
g2 <- g2 + theme_classic() +
  labs(title = "", 
       x = "", 
       y = "Mean Intensity")

print(g2)

#-------

g3 <- ggplot(detected_counts, aes(x = Tissue, y = DetectedCount, fill = Turbo)) 
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

# convert human gene names to mouse

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

# data <- data %>%
#   separate(genes, into = paste0("genes"), sep = ";", remove = TRUE) 

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

match_index <- match(colnames(scaled_expression), metadata$sampleID)
metadata_reordered <- metadata[match_index, ]

tissue_list <- as.factor(paste(metadata_reordered$Tissue,"-",metadata_reordered$Turbo))

# Create a color mapping for metadata

tissue_colors <- c("DRG - T" = "#f55c3a",
                   "DRG - TC" = "#f1ccc4",
                   "LSC - T" = "#edb127",
                   "LSC - TC" = "#f7e1ae",
                   "SCN - T" = "#29c99d",
                   "SCN - TC" = "#abcfc5",
                   "paw - T" = "#95459b",
                   "paw - TC" = "#e7d9e8"
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

PATH_results = "./output/"

pdf(file = paste0(PATH_results, "/neuronal-heatmap.pdf"), height = 5, width = 8)
draw(ht_list, heatmap_legend_side = "right")
dev.off()

#-------------------------------------------------------------------------------

# export data for downstream analyses

data <- df
rownames(data) <- data$GeneID
data$GeneID <- NULL

colData <- metadata

data_trans  <- as.data.frame(log2(data))
ranked_data <- data_trans %>% mutate_all(rank)

merged_data <- data_trans %>%
  tibble::rownames_to_column(var = "proteins") %>%
  pivot_longer(cols = -proteins, names_to = "sampleID", values_to = "Log2Intensity") %>%
  left_join(ranked_data %>% 
              tibble::rownames_to_column(var = "proteins") %>%
              pivot_longer(cols = -proteins, names_to = "sampleID", values_to = "Rank"),
            by = c("proteins", "sampleID")) %>%
  left_join(colData, by = "sampleID") %>%
  separate(proteins, into = "genes", sep = ";", remove = FALSE) %>%
  mutate(genes = trimws(as.character(genes)))

df <- merged_data %>%
  group_by(proteins, genes, sampleID) %>%
  summarize(expression = mean(Log2Intensity, na.rm = TRUE))

df <- df %>% pivot_wider(names_from = sampleID, values_from = expression)
df <- as.data.frame(df)

df <- df %>%
  mutate(proteins = gsub("\\n", "", proteins)) %>%
  distinct(proteins, .keep_all = TRUE)

head(df)

# reduce colData for merge
df_meta <- colData
df_meta <- df_meta[!duplicated(df_meta), ]

write.csv(df,      "./data/matrix-for-limma.csv", na="NA", eol = "\n", row.names = FALSE)
write.csv(df_meta, "./data/colData-for-limma.csv", na="NA", eol = "\n", row.names = FALSE)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Plot Dynamic Range fpr each Tissue + Turbo Condition

ranked_mean <- merged_data %>%
  select(proteins, genes, Tissue, Turbo, Log2Intensity) %>%
  group_by(proteins, Tissue, Turbo, genes) %>%
  summarize(mean_Log2Intensity = mean(Log2Intensity, na.rm = TRUE), .groups = 'drop') %>%
  group_by(Tissue, Turbo) %>%
  mutate(Rank = rank(-mean_Log2Intensity, ties.method = "first")) %>%
  ungroup()

# superimpose key genes on plot
named_mean <- ranked_mean %>%
  filter(genes %in% neurons) %>%
  group_by(genes, Tissue, Turbo) %>%
  filter(Rank == min(Rank)) %>%
  ungroup()

# head(named_mean)

#-------

ranked_mean <- na.omit(ranked_mean)
named_mean  <- na.omit(named_mean)

# dynamic range by group
g <- ggplot() +
  geom_line(data = ranked_mean, mapping = aes(x = Rank, y = mean_Log2Intensity), color = "#FFE6BF", linewidth = 2) +
  geom_point(data = named_mean, mapping = aes(x = Rank, y = mean_Log2Intensity), size = 2, shape = 21, fill = "#21130d", colour = "#21130d")
g <- g + ggrepel::geom_label_repel(data = named_mean, mapping = aes(x = Rank, y = mean_Log2Intensity, label = genes), 
                                   stat = "identity", color = "#21130d", segment.color = 'grey50', 
                                   force = 50, box.padding = 0.35, point.padding = 0.5, max.overlaps = 45) 
g <- g + facet_grid(Tissue ~ Turbo) + theme_bw() 
g <- g + labs(title = "Key proteins", x = "Rank", y = "Log2 Intensity")

print(g)

pdf(file = paste(PATH_results, "dynamic-range.pdf", sep=""), width = 9, height = 6)
print(g)
dev.off()

#-------

#plot only TurboID positive samples
ranked_mean <- ranked_mean[ranked_mean$Turbo %in% "T", ]
named_mean  <- named_mean[named_mean$Turbo %in% "T", ]

g <- ggplot() +
  geom_line(data = ranked_mean, mapping = aes(x = Rank, y = mean_Log2Intensity), color = "#FFE6BF", linewidth = 2) +
  geom_point(data = named_mean, mapping = aes(x = Rank, y = mean_Log2Intensity), size = 2, shape = 21, fill = "#21130d", colour = "#21130d")
g <- g + ggrepel::geom_label_repel(data = named_mean, mapping = aes(x = Rank, y = mean_Log2Intensity, label = genes), 
                                   stat = "identity", color = "#21130d", segment.color = 'grey50', 
                                   force = 50, box.padding = 0.35, point.padding = 0.5, max.overlaps = 45) 
g <- g + facet_grid(Tissue ~ .) + theme_bw() 
g <- g + labs(title = "Key proteins", x = "Rank", y = "Log2 Intensity")

print(g)

pdf(file = paste(PATH_results, "dynamic-range_TurboID.pdf", sep=""), width = 6, height = 7)
print(g)
dev.off()

#-------------------------------------------------------------------------------

df      <- read.csv("./data/matrix-for-limma.csv", header = TRUE)
colData <- read.csv("./data/colData-for-limma.csv", header = TRUE)

head(df)


