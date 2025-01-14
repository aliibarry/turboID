library(dplyr)
library(tidyr)
library(limma)
library(ggplot2)
library(VennDiagram)

df      <- read.csv("./data/matrix-for-limma.csv", header = TRUE)
colData <- read.csv("./data/colData-for-limma.csv", header = TRUE)

PATH_results = "./output/"

#-------------------------------------------------------------------------------

tissues <- unique(colData$Tissue)
tissue_counts <- list()

# iterate over tissues to count overlapping PGs for T/TC
for (tissue in tissues) {
  tissue_samples <- colData[colData$Tissue == tissue, ]
  turbo_groups <- split(tissue_samples$sampleID, tissue_samples$Turbo)
  tissue_data <- df[, tissue_samples$sampleID, drop = FALSE]
  
  presence_T  <- !is.na(tissue_data[, turbo_groups$T, drop = FALSE])
  presence_TC <- !is.na(tissue_data[, turbo_groups$TC, drop = FALSE])
  
  count_T       <- sum(rowSums(presence_T) > 0 & rowSums(presence_TC) == 0)  # Specific to T
  count_TC      <- sum(rowSums(presence_TC) > 0 & rowSums(presence_T) == 0) # Specific to TC
  count_overlap <- sum(rowSums(presence_T) > 0 & rowSums(presence_TC) > 0) # Overlapping
  
  tissue_counts[[tissue]] <- data.frame(
    Category = c("T", "TC", "Overlap"),
    Count = c(count_T, count_TC, count_overlap)
  )
}

plot_data <- do.call(rbind, lapply(names(tissue_counts), function(tissue) {
  cbind(Tissue = tissue, tissue_counts[[tissue]])
}))

g <- ggplot(plot_data, aes(x = Tissue, y = Count, fill = Category)) 
g <- g + geom_bar(stat = "identity", position = "stack") 
g <- g + geom_text(
  aes(label = Count),
  position = position_stack(vjust = 0.5),
  color = "black") 
g <- g + labs(title = "", x = "", y = "PG Count") 
g <- g + theme_bw() +
  scale_fill_manual(
    values = c("T" = "#afc9e3", "TC" = "#518ac4", "Overlap" = "#74a9cf"),
    name = "Category")

print(g)

pdf(file = paste(PATH_results, "PG_counts_by_tissue.pdf", sep=""), width = 5, height = 4)
print(g)
dev.off()

#-------------------------------------------------------------------------------

#dir.create("./output/enrichments/lm_all")
PATH_results = "./output/enrichments/lm_all/"

mat <- df %>% distinct(genes, .keep_all = TRUE) #6109 distinct genes
rownames(mat) <- mat$genes

mat$proteins <- NULL
mat$genes    <- NULL

# reorder colData to match matrix
colData <- colData[colData$sampleID %in% colnames(mat), ]

index   <- match(colnames(mat), colData$sampleID)
colData <- colData[index, ]

# filter for 90% present (ie. one exerimental group)
# mat[is.na(mat)] <- 0
NAs <- rowMeans(is.na(mat))
mat <- mat[NAs <= 0.90, ] #0.09756098 (8/82 samples == smallest category)
dim(mat)

# colData MUST match order of mat
head(colnames(mat))
head(colData$sampleID)

colData$TissueTurbo <- paste0(colData$Tissue, "_", colData$Turbo)

design <- model.matrix(~ 0 + TissueTurbo, data = colData) #select design
head(design)

fit <- lmFit(mat, design)

colnames(fit) #check possible comparisons

# Run hypothesis testing, adjust tissue as needed
contrast_matrix <- makeContrasts( TissueTurboDRG_T -  TissueTurboDRG_TC, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust="BH", number=Inf)

head(results)

sig_de <- results[results$adj.P.Val < 0.05, , drop = FALSE] #extract significant proteins if there
head(sig_de)

sig_de <- na.omit(sig_de)

# volcano plotting
res       <- as.data.frame(results) 
mutateddf <- mutate(res, Sig=ifelse(res$adj.P.Val<0.05 & abs(res$logFC)>1, "FDR < 0.05 & LFC > 1", ifelse("NS")))
input     <- cbind(gene=rownames(res), mutateddf) 

volc = ggplot(input, aes(logFC, -log10(P.Value))) + geom_point(aes(col=Sig)) +
  #scale_color_manual(values = c("#0D0887FF","#9512A1FF", "grey")) + 
  scale_color_manual(values = c("#B63679ff", "grey")) + 
  ggrepel::geom_text_repel(data=subset(input, input$gene %in% rownames(sig_de)),
                           aes(label=gene), size=4, segment.alpha= 0.2, force =2, max.overlaps=16) 
volc <- volc + theme_bw() + theme(aspect.ratio=1)
volc <- volc + theme(legend.position="bottom", axis.text.y = element_text(size= 12, face="bold"), 
                     axis.title.y = element_text(size=14), axis.title.x = element_text(size= 14), 
                     axis.text.x = element_text(size= 12), legend.title=element_text(size=14), 
                     legend.text=element_text(size=14), plot.title=element_text(size=12, hjust = 0.5)) + 
  ggtitle("DRG")

print(volc)

pdf(file = paste0(PATH_results, "volcano-DRG.pdf"), height = 4, width = 4)
print(volc)
dev.off()

pdf(file = paste0(PATH_results, "volcano-DRG_big.pdf"), height = 6, width = 6)
print(volc)
dev.off()

table(sig_de$logFC > 1)
write.csv(results, paste(PATH_results, "DEP-analysis-limma_DRG.csv"))

#-------------------------------------------------------------------------------

# load&filter limma results above above run on each tissue
drg <- read.csv("./output/enrichments/lm_all/ DEP-analysis-limma_DRG.csv", header = TRUE, row.names = 1)
scn <- read.csv("./output/enrichments/lm_all/ DEP-analysis-limma_SCN.csv", header = TRUE, row.names = 1)
lsc <- read.csv("./output/enrichments/lm_all/ DEP-analysis-limma_LSC.csv", header = TRUE, row.names = 1)
paw <- read.csv("./output/enrichments/lm_all/ DEP-analysis-limma_paw.csv", header = TRUE, row.names = 1)

drg <- na.omit(drg)
scn <- na.omit(scn)
lsc <- na.omit(lsc)
paw <- na.omit(paw)

head(drg)
str(colData)

# extract "significant" enrichments where gene present in both Turbo and control conditions (LFC > 1, FDR < 0.05)
drg_de <- data.frame(Gene = rownames(drg)[drg$logFC > 1 & drg$adj.P.Val < 0.05], Tissue = "DRG", Type = "DEP")
scn_de <- data.frame(Gene = rownames(scn)[scn$logFC > 1 & scn$adj.P.Val < 0.05], Tissue = "SCN", Type = "DEP")
lsc_de <- data.frame(Gene = rownames(lsc)[lsc$logFC > 1 & lsc$adj.P.Val < 0.05], Tissue = "LSC", Type = "DEP")
paw_de<- data.frame(Gene = rownames(paw)[paw$logFC > 1 & paw$adj.P.Val < 0.05], Tissue = "paw", Type = "DEP")

all_degs <- bind_rows(drg_de, scn_de, lsc_de, paw_de)
str(all_degs)

#-------------------------------------------------------------------------------

# Calculate enrichments
# NOTE, original calculations from Julia on PGS (not disinct genes), and then filtered after (ie. slightly elevated pre-filter)

PATH_results = "./output/enrichments/"

df      <- read.csv("./data/matrix-for-limma.csv", header = TRUE)
colData <- read.csv("./data/colData-for-limma.csv", header = TRUE)

df <- df %>% distinct(genes, .keep_all = TRUE)

rownames(df) <- df$genes
df$proteins  <- NULL
df$genes     <- NULL

# confirm matching order
colData <- colData[colData$sampleID %in% colnames(df), ]
index   <- match(colnames(df), colData$sampleID)
colData <- colData[index, ]

# filter genes by 75% TurboID, 0% Turbo Control
filter_genes <- function(tissue, colData, df) {
  tissue_data <- colData %>% filter(Tissue == tissue)

  t_samples  <- df[, tissue_data$sampleID[tissue_data$Turbo == "T"], drop = FALSE]
  tc_samples <- df[, tissue_data$sampleID[tissue_data$Turbo == "TC"], drop = FALSE]
  
  genes_to_keep <- rownames(df)[
      rowMeans(!is.na(t_samples)) >= 0.75 &  # At least 75% of T samples are non-NA
      rowSums(!is.na(tc_samples), na.rm = TRUE) == 0  # All TC samples must be NA
  ]
  
  data.frame(Tissue = tissue, Gene = genes_to_keep)
}

# Apply the function for each tissue
drg_genes <- filter_genes("DRG", colData, df)
lsc_genes <- filter_genes("LSC", colData, df)
scn_genes <- filter_genes("SCN", colData, df)
paw_genes <- filter_genes("paw", colData, df)

# Combine results into a single data frame
all_filtered      <- bind_rows(drg_genes, lsc_genes, scn_genes, paw_genes)
all_filtered$Type <- "75% Filter"

# Merge DEGs and filtered genes into one data frame
merged_genes <- bind_rows(all_degs, all_filtered)

# View the merged result
head(merged_genes)

#-------------------------------------------------------------------------------

# Venn Plotting
drg_genes <- unique(merged_genes$Gene[merged_genes$Tissue == "DRG"])
scn_genes <- unique(merged_genes$Gene[merged_genes$Tissue == "SCN"])
lsc_genes <- unique(merged_genes$Gene[merged_genes$Tissue == "LSC"])
paw_genes <- unique(merged_genes$Gene[merged_genes$Tissue == "paw"])

venn.plot <- venn.diagram(
  x = list(DRG = drg_genes, SCN = scn_genes, LSC = lsc_genes, paw = paw_genes),
  category.names = c("DRG", "SCN", "LSC", "paw"),
  filename = NULL,
  output = TRUE,
  imagetype = "pdf",
  height = 300, width = 300, resolution = 300,
  col = "black",
  fill = c("#f55c3a", "#29c99d", "#edb127", "#95459b"),
  alpha = 0.3,
  label.col = "black",
  cex = 1,
  fontface = "plain",
  # main.fontfamily = "plain",
  #sub.fontface = "normal",
  #fontfamily = "plain",
  cat.cex = 1
)

grid.newpage()
grid::grid.draw(venn.plot)

pdf(file = paste0(PATH_results, "venn.pdf"), height = 6, width = 6)
grid::grid.draw(venn.plot)
dev.off()

# Box plotting
plot_data <- merged_genes %>%
  dplyr::group_by(Tissue, Type) %>%
  dplyr::summarise(Count = n()) %>%
  dplyr::ungroup()

g <- ggplot(plot_data, aes(x = Tissue, y = Count, fill = Type)) 
g <- g + geom_bar(stat = "identity", position = "stack") 
g <- g + geom_text(aes(label = Count),
  position = position_stack(vjust = 0.5),
  color = "black") 
g <- g + labs(title = "Turbo Enrichments", x = " ", y = "Gene Count") 
g <- g + theme_bw() +
  scale_fill_manual(
    values = c("DEP" = "#afc9e3", "75% Filter" = "#518ac4"),
    name = "Type")

print(g)

pdf(file = paste0(PATH_results, "enrichments.pdf"), width = 5, height = 4)
print(g)
dev.off()

write.csv(merged_genes, "./output/enrichments.csv")
