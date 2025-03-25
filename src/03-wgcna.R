library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(ComplexHeatmap)

library(WGCNA)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)

#dir.create("./output/WGCNA/")
PATH_results = "./output/WGCNA/"

df          <- read.csv("./data/matrix-for-limma.csv", header = TRUE)
colData     <- read.csv("./data/colData-for-limma.csv", header = TRUE)
enrichments <- read.csv("./output/enrichments_75filt.csv", 
                        check.names = FALSE, header = TRUE, row.names = 1)
#-------------------------------------------------------------------------------                                 

# test_df <- enrichments$Gene[enrichments$Tissue == "paw"] 
test_df <- enrichments$Gene[!duplicated(enrichments$Gene)]

head(test_df)

df_filt <- df[df$genes %in% test_df, ]
df_filt <- df_filt[!duplicated(df_filt$genes), ]

rownames(df_filt) <- df_filt$genes

df_filt$genes    <- NULL
df_filt$proteins <- NULL

df_filt <- df_filt[colnames(df_filt) %in% colData$sampleID[colData$Turbo == "T"]]
colData <- colData[match(colnames(df_filt), colData$sampleID), ]
rownames(colData) <- colData$sampleID

# #keep genes in more than one tissue only
sample_tissue_map <- setNames(colData$Tissue, colData$sampleID)
sample_tissues <- sample_tissue_map[colnames(df_filt)]

filter_genes <- function(expression_values) {
  expressed_tissues <- unique(na.omit(sample_tissues[!is.na(expression_values)]))
  length(expressed_tissues) > 1
}

# Apply the filter function to rows (genes)
df_filt <- df_filt[apply(df_filt, 1, filter_genes), ]
 
# min_value    <- min(df_filt, na.rm = TRUE)
# median_value <- median(as.matrix(df_filt), na.rm = TRUE)
# 
# impute_median_per_group <- function(df, colData, group_col) {
# 
#   colData   <- colData[colnames(df), , drop = FALSE]
#   group_ids <- colData[[group_col]]
# 
#   # Apply function to each row (gene)
#   df <- t(apply(df, 1, function(gene_expr) {
# 
#     group_medians <- tapply(gene_expr, group_ids, median, na.rm = TRUE)
#     overall_median <- median(gene_expr, na.rm = TRUE)
#     gene_expr[is.na(gene_expr)] <- ifelse(is.na(group_medians[group_ids[is.na(gene_expr)]]),
#                                           overall_median, group_medians[group_ids[is.na(gene_expr)]])
#     return(gene_expr)
#   }))
# 
#   return(df)
# }
# 
# df_filt <- impute_median_per_group(df_filt, colData, "Tissue")

df_filt <- as.data.frame(df_filt)

table(is.na(df_filt))

#filter for more variable genes
gene_variances <- apply(df_filt, 1, var, na.rm = TRUE)
hist(gene_variances)
range(gene_variances)

df_filt <- df_filt[gene_variances > 0, ]
dim(df_filt)

#-------------------------------------------------------------------------------

# options(stringsAsFactors = FALSE)
# enableWGCNAThreads(nThreads = 8)

datExpr <- as.data.frame(t(df_filt))

powers <- c(1:50)  # Range of powers to test
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed")

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], type="o",
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit", 
     main="Soft Threshold Selection")
abline(h=0.8, col="red")  # Choose power where fit is ~0.9

power <- 14  # Adjust based on previous step
adjacency <- adjacency(datExpr, power = power, type = "signed")

# Convert adjacency into a Topological Overlap Matrix (TOM)
TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1 - TOM

#-------------------------------------------------------------------------------

set.seed(52)
selectedGenes <- sample(ncol(datExpr), 400)  
TOMsubset <- TOM[selectedGenes, selectedGenes]

# Visualize with a heatmap
heatmap(TOMsubset, col = viridis::viridis(100), symm = TRUE)

pheatmap(TOMsubset, 
         color = viridis::viridis(100), 
         clustering_method = "average", 
         main = "TOM Heatmap (Subset of Genes)")

#-------------------------------------------------------------------------------


# Hierarchical clustering of genes
geneTree <- hclust(as.dist(dissTOM), method = "average")

# Plot the dendrogram
plot(geneTree, main = "Gene Clustering Dendrogram", sub = "", xlab = "")

pdf(file = paste(PATH_results, "dendrogram_raw.pdf", sep=""), width = 8, height = 8)
plot(geneTree, main = "Gene Clustering Dendrogram", sub = "", xlab = "")
dev.off()

# Dynamic tree cut to identify modules
minModuleSize <- 30  # Minimum module size
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)

# Convert module labels to colors
moduleColors <- labels2colors(dynamicMods)
table(moduleColors)

# Plot dendrogram with module colors
plotDendroAndColors(geneTree, moduleColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

pdf(file = paste(PATH_results, "dendrogram_colours.pdf", sep=""), width = 8, height = 5)
plotDendroAndColors(geneTree, moduleColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03, addGuide = FALSE, guideHang = 0.05)

dev.off()

#-------------------------

traits <- colData
rownames(traits) <- traits$sampleID

traits$Turbo  <- NULL
traits$Sample <- NULL
traits$ID     <- NULL
traits$sampleID <- NULL
traits$Tissue <- as.integer(as.factor(traits$Tissue))
traits$Sex    <- as.integer(as.factor(traits$Sex))

# Check and match samples
traits <- traits[match(rownames(datExpr), rownames(traits)),]

# Correlate module eigengenes with traits
MEs <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
MEs <- MEs[, !colnames(MEs) %in% "MEgrey"]

moduleTraitCor <- cor(MEs, traits, use = "p", method = "spearman")
moduleTraitP   <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# extrat significant correlations
p_values <- moduleTraitP[, 1]
adjusted_p_values <- p.adjust(p_values, method = "BH")
significant_names <- rownames(moduleTraitP)[adjusted_p_values < 0.05]

significant_names

# Plot heatmap of module-trait relationships
labeledHeatmap(Matrix = moduleTraitCor, xLabels = colnames(traits),
               yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE,
               textMatrix = signif(moduleTraitCor, 2),
               main = "Module-Trait Relationships")

tissue_values <- moduleTraitCor[, "Tissue"]
ordered_rows  <- order(tissue_values, decreasing = TRUE)
moduleTraitCor_ordered <- moduleTraitCor[ordered_rows, , drop = FALSE]
moduleTraitP_ordered    <- moduleTraitP[ordered_rows, , drop = FALSE]

pdf(file = paste(PATH_results, "modules-heatmap.pdf", sep=""), width = 6, height = 7)
labeledHeatmap(Matrix = moduleTraitCor_ordered, 
               xLabels = colnames(traits),
               yLabels = rownames(moduleTraitCor_ordered), 
               ySymbols = rownames(moduleTraitCor_ordered), 
               colorLabels = FALSE,
               textMatrix = signif(moduleTraitCor_ordered, 2),
               main = "Module-Trait Relationships",
               colors = viridis(100))  # Use the viridis color scale with 100 shades
dev.off()

str(moduleTraitP)

#---------------------------

traits
tissue_names <- c("DRG","LSC", "paw", "SCN") #adjust order as needed, check names
traits$Tissue_name <- tissue_names[traits$Tissue] 

MEs_matrix <- t(MEs)  # Transpose MEs, so rows are modules and columns are samples

match_index <- match(colnames(MEs_matrix), rownames(traits))
traits <- traits[match_index, ]

rownames(MEs_matrix) <- names(MEs)  # Rows are module names

tail(MEs_matrix)
tail(traits)

MEs_matrix_sig <- MEs_matrix[rownames(MEs_matrix) %in% significant_names,]

match_index <- match(colnames(MEs_matrix_sig), rownames(traits))
traits <- traits[match_index, ]

tissue_list <- as.factor(traits$Tissue_name)

# Create a color mapping for metadata
tissue_colors <- c("DRG" = "#f55c3a",
                   "LSC" = "#edb127",
                   "SCN" = "#29c99d",
                   "paw" = "#95459b"
)

# Assign colors to metadata levels
col_fun <- tissue_colors[tissue_list]

Heatmap(MEs_matrix_sig, 
        name = "Module Eigengenes", 
        column_title = "Tissue",
        row_title = "Module",
        show_row_names = TRUE, 
        show_column_names = TRUE,  # Show tissue names on top
        cluster_rows = TRUE,  # Cluster the rows (modules)
        cluster_columns = TRUE,  # Cluster the columns (samples)
        top_annotation = HeatmapAnnotation(tissue = tissue_list, col = list(tissue = col_fun)),
        heatmap_legend_param = list(title = "Eigengene Value")
)

pdf(file = paste(PATH_results, "modules-by-tissue-heatmap.pdf", sep=""), width = 7, height = 6)
Heatmap(MEs_matrix_sig, 
        name = "Module Eigengenes", 
        column_title = "Tissue",
        row_title = "Module",
        show_row_names = TRUE, 
        show_column_names = FALSE,  # Show tissue names on top
        cluster_rows = TRUE,  # Cluster the rows (modules)
        cluster_columns = TRUE,  # Cluster the columns (samples)
        top_annotation = HeatmapAnnotation(tissue = tissue_list, col = list(tissue = col_fun)),
        heatmap_legend_param = list(title = "Eigengene Value")
)

dev.off()

#-------------------------------------------------------------------------------

#Check MEs by tissue
MEs_long <- as.data.frame(MEs_matrix_sig) %>%
  mutate(Module = rownames(MEs_matrix_sig)) %>%
  pivot_longer(cols = -Module, names_to = "Sample", values_to = "ME_value")

MEs_long <- MEs_long %>%
  left_join(colData, by = c("Sample" = "sampleID"))

MEs_long$Tissue <- factor(MEs_long$Tissue, levels = c("paw", "SCN", "DRG", "LSC"))

g <- ggplot(MEs_long, aes(x = Tissue, y = ME_value, fill = Tissue)) 
g <- g + geom_boxplot() 
g <- g + theme_bw() + facet_wrap(~ Module, scales = "fixed") 
g <- g +  # Facet by module, with independent y-scales
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  labs(x = " ", y = "ME Value", title = "MEs by Tissue and Module")

print(g)

pdf(file = paste(PATH_results, "modules-by-tissue.pdf", sep=""), width = 8, height = 7)
print(g)
dev.off()

#-------------------------------------------------------------------------------

module <- "green"
moduleGenes <- colnames(datExpr)[moduleColors == module]

# Calculate module membership (MM) scores
MM <- cor(datExpr[, moduleGenes], MEs[, paste0("ME",module)], use = "p")
names(MM) <- row.names(MM)

topHubGenes <- names(sort(MM, decreasing = TRUE)[1:10])
print(topHubGenes)

MM <- as.data.frame(MM)

write.csv(MM, paste0(PATH_results,"module_green.csv"))

#-------------------------------------------------------------------------------

# run GO analysis on each module
background <- test_df
# background <- rownames(df_filt)

module <- "green"
moduleGenes <- colnames(datExpr)[moduleColors == module]

rm(ego2, ego)
ego <- enrichGO(gene          = moduleGenes,
                universe      = background,
                OrgDb         = org.Mm.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

# remove redundancy in the GO terms
ego2    <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min, measure = 'Wang')

pdf(file = paste(PATH_results, "BP_green-network.pdf", sep=""), width = 8, height = 8)
goplot(ego2)
dev.off()

pdf(file = paste(PATH_results, "BP_green-barplot.pdf", sep=""), width = 6, height = 3)
mutate(ego2, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")
dev.off()

pdf(file = paste(PATH_results, "BP_green-upset.pdf", sep=""), width = 7, height = 4)
upsetplot(ego2)
dev.off()

write.csv(ego2, paste(PATH_results, "BP_green.csv"))

#-----

ego <- enrichGO(gene          = moduleGenes,
                universe      = background,
                OrgDb         = org.Mm.eg.db,
                keyType       = "SYMBOL",
                ont           = "CC", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

# remove redundancy in the GO terms
ego2    <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min, measure = 'Wang')

pdf(file = paste(PATH_results, "CC_green-network.pdf", sep=""), width = 8, height = 8)
goplot(ego2)
dev.off()

pdf(file = paste(PATH_results, "CC_green-barplot.pdf", sep=""), width = 6, height = 3)
mutate(ego2, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")
dev.off()

pdf(file = paste(PATH_results, "CC_green-upset.pdf", sep=""), width = 7, height = 4)
upsetplot(ego2)
dev.off()

write.csv(ego2, paste(PATH_results, "CC_green.csv"))

#---

ego <- enrichGO(gene          = moduleGenes,
                universe      = background,
                OrgDb         = org.Mm.eg.db,
                keyType       = "SYMBOL",
                ont           = "MF", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

# # remove redundancy in the GO terms
# ego2    <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min, measure = 'Wang')

pdf(file = paste(PATH_results, "MF_green-network.pdf", sep=""), width = 8, height = 8)
goplot(ego)
dev.off()

pdf(file = paste(PATH_results, "MF_green-barplot.pdf", sep=""), width = 6, height = 3)
mutate(ego, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")
dev.off()

pdf(file = paste(PATH_results, "MF_green-upset.pdf", sep=""), width = 7, height = 4)
upsetplot(ego)
dev.off()

write.csv(ego, paste(PATH_results, "MF_green.csv"))

#---



print(g)

#-------------------------------------------------------------------------------

#SCRATCH



