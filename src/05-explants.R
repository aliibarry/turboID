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

dir.create("./output/explants/QC")
PATH_results = "./output/explants/QC/"

df      <- read.csv("./data/explants/combined_genes.csv", row.names = 1, check.names = FALSE)
colData <- read.csv("./data/explants/colData.csv", row.names = 1)

#-------------------------------------------------------------------------------

df$GeneID <- rownames(df)

df_long <- df %>%
  pivot_longer(cols = -GeneID, names_to = "sampleID", values_to = "Expression") %>%
  left_join(colData, by = "sampleID")

summed_expression <- df_long %>%
  group_by(sampleID, Condition, Turbo) %>%
  summarize(SumCount = sum(Expression, na.rm = TRUE), .groups = "drop") %>%
  mutate(Condition = factor(Condition, levels = c("Ox", "V")))

mean_intensity <- df_long %>%
  group_by(sampleID, Condition, Turbo) %>%
  summarize(MeanIntensity = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
  mutate(Condition = factor(Condition, levels = c("Ox", "V")))

detected_counts <- df_long %>%
  group_by(sampleID, Condition, Turbo) %>%
  summarize(DetectedCount = sum(!is.na(Expression)), .groups = "drop") %>%
  mutate(Condition = factor(Condition, levels = c("Ox", "V")))

#-------

g1 <- ggplot(summed_expression, aes(x = Condition, y = SumCount, fill = Turbo)) 
g1 <- g1 + geom_boxplot() 
g1 <- g1 + theme_classic() +
  labs(title = "", 
       x = "", 
       y = "Summed Expression")

print(g1)

#------

g2 <- ggplot(mean_intensity, aes(x = Condition, y = MeanIntensity, fill = Turbo)) 
g2 <- g2 + geom_boxplot() 
g2 <- g2 + theme_classic() +
  labs(title = "", 
       x = "", 
       y = "Mean Intensity")

print(g2)

#-------

g3 <- ggplot(detected_counts, aes(x = Condition, y = DetectedCount, fill = Turbo)) 
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

sex    <- as.factor(colData$Sex)
tissue <- as.factor(colData$Condition)
turbo  <- as.factor(colData$Turbo)

g1 <- ggbiplot(pca, choices = c(1,2), scale = 0.5,
               groups = interaction(tissue),
               #ellipse = TRUE,
               ellipse.prob = 0.95,
               labels = NULL,
               point.size = 4,
               labels.size = 4, alpha = 1, var.axes = FALSE,
               circle  = TRUE, circle.prob = 0.5,
               varname.size = 3,
               varname.adjust = 1.5,
               varname.abbrev = FALSE)

g1 <- g1 + theme(legend.position = 'right') + theme_bw()

g2 <- ggbiplot(pca, choices = c(1,2), scale = 0.5,
               groups = interaction(turbo, tissue),
               #ellipse = TRUE,
               ellipse.prob = 0.95,
               point.size = 4,
               labels = NULL,
               labels.size = 4, alpha = 1, var.axes = FALSE,
               circle  = TRUE, circle.prob = 0.5,
               varname.size = 3,
               varname.adjust = 1.5,
               varname.abbrev = FALSE)

g2 <- g2 + theme(legend.position = 'right') + theme_bw()

print(g2)
print(g1)

g3 <- ggbiplot(pca, choices = c(1,2), scale = 0.5,
               groups = interaction(turbo),
               #ellipse = TRUE,
               ellipse.prob = 0.95,
               point.size = 4,
               labels = NULL,
               labels.size = 4, alpha = 1, var.axes = FALSE,
               circle  = TRUE, circle.prob = 0.5,
               varname.size = 3,
               varname.adjust = 1.5,
               varname.abbrev = FALSE)

g3 <- g3 + theme(legend.position = 'right') + theme_bw()

g4 <- ggbiplot(pca, choices = c(1,2), scale = 0.5,
               groups = interaction(sex),
               #ellipse = TRUE,
               ellipse.prob = 0.95,
               point.size = 4,
               labels = NULL,
               labels.size = 4, alpha = 1, var.axes = FALSE,
               circle  = TRUE, circle.prob = 0.5,
               varname.size = 3,
               varname.adjust = 1.5,
               varname.abbrev = FALSE)

g4 <- g4 + theme(legend.position = 'right') + theme_bw()

PATH_results = "./output/explants/"

pdf(file = paste(PATH_results, "pca_all.pdf", sep=""), width = 8, height = 6)
grid.arrange(g1, g2, g3, g4, ncol = 2)
dev.off()

#-------------------------------------------------------------------------------

# remove TurboID controls first

colData <- colData[colData$Turbo %in% "T", ]

mat <- as.matrix(df[, which(colnames(df) %in% colData$sampleID)])
mat <- mat[complete.cases(mat), ]

pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)

index   <- match(colnames(mat), colData$sampleID)
colData <- colData[index, ]

sex    <- as.factor(colData$Sex)
tissue <- as.factor(colData$Condition)
turbo  <- as.factor(colData$Turbo)

g1 <- ggbiplot(pca, choices = c(1,2), scale = 0.5,
               groups = interaction(tissue),
               #ellipse = TRUE,
               ellipse.prob = 0.95,
               labels = NULL,
               point.size = 4,
               labels.size = 4, alpha = 1, var.axes = FALSE,
               circle  = TRUE, circle.prob = 0.5,
               varname.size = 3,
               varname.adjust = 1.5,
               varname.abbrev = FALSE)

g1 <- g1 + theme(legend.position = 'right') + theme_bw()
#g1 <- g1 + ggtitle("Sample subset")

g2 <- ggbiplot(pca, choices = c(1,2), scale = 0.5, 
               groups = interaction(turbo),
               #ellipse = TRUE,
               ellipse.prob = 0.95,
               point.size = 4,
               labels = NULL,
               labels.size = 4, alpha = 1, var.axes = FALSE,
               circle  = TRUE, circle.prob = 0.5,
               varname.size = 3,
               varname.adjust = 1.5,
               varname.abbrev = FALSE)

g2 <- g2 + theme(legend.position = 'right') + theme_bw()
#g2 <- g2 + ggtitle("Cohort")
print(g2)
print(g1)

g3 <- ggbiplot(pca, choices = c(1,2), scale = 0.5,
               groups = interaction(sex, tissue),
               #ellipse = TRUE,
               ellipse.prob = 0.95,
               point.size = 4,
               labels = NULL,
               labels.size = 4, alpha = 1, var.axes = FALSE,
               circle  = TRUE, circle.prob = 0.5,
               varname.size = 3,
               varname.adjust = 1.5,
               varname.abbrev = FALSE)

g3 <- g3 + theme(legend.position = 'right') + theme_bw()

g4 <- ggbiplot(pca, choices = c(1,2), scale = 0.5,
               groups = interaction(sex),
               #ellipse = TRUE,
               ellipse.prob = 0.95,
               point.size = 4,
               labels = NULL,
               labels.size = 4, alpha = 1, var.axes = FALSE,
               circle  = TRUE, circle.prob = 0.5,
               varname.size = 3,
               varname.adjust = 1.5,
               varname.abbrev = FALSE)

g4 <- g4 + theme(legend.position = 'right') + theme_bw()


pdf(file = paste(PATH_results, "pca.pdf", sep=""), width = 8, height = 6)
grid.arrange(g1, g2, g3, g4, ncol = 2)
dev.off()

