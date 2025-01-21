library(dplyr)
library(tidyr)
library(limma)
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

#dir.create("./output/explants/QC")
PATH_results = "./output/explants/"

df      <- read.csv("./data/explants/combined_genes.csv", row.names = 1, check.names = FALSE)
colData <- read.csv("./data/explants/colData.csv", row.names = 1)

#-------------------------------------------------------------------------------

head(df[1:5])

head(colData)

#calculate enrichments as before
#compare to enrichment lists (DRG-specific & combined)

Conditions <- unique(colData$Condition)
Condition_counts <- list()

# iterate over Conditions to count overlapping PGs for T/TC
for (Condition in Conditions) {
  Condition_samples <- colData[colData$Condition == Condition, ]
  turbo_groups <- split(Condition_samples$sampleID, Condition_samples$Turbo)
  Condition_data <- df[, Condition_samples$sampleID, drop = FALSE]
  
  presence_T  <- !is.na(Condition_data[, turbo_groups$T, drop = FALSE])
  presence_TC <- !is.na(Condition_data[, turbo_groups$TC, drop = FALSE])
  
  count_T       <- sum(rowSums(presence_T) > 0 & rowSums(presence_TC) == 0)  # Specific to T
  count_TC      <- sum(rowSums(presence_TC) > 0 & rowSums(presence_T) == 0) # Specific to TC
  count_overlap <- sum(rowSums(presence_T) > 0 & rowSums(presence_TC) > 0) # Overlapping
  
  Condition_counts[[Condition]] <- data.frame(
    Category = c("T", "TC", "Overlap"),
    Count = c(count_T, count_TC, count_overlap)
  )
}

plot_data <- do.call(rbind, lapply(names(Condition_counts), function(Condition) {
  cbind(Condition = Condition, Condition_counts[[Condition]])
}))

g <- ggplot(plot_data, aes(x = Condition, y = Count, fill = Category)) 
g <- g + geom_bar(stat = "identity", position = "stack") 
g <- g + geom_text(
  aes(label = Count),
  position = position_stack(vjust = 0.5),
  color = "black") 
g <- g + labs(title = "", x = "", y = "PG Count") 
g <- g + theme_bw() +
  scale_fill_manual(
    values = c("T" = "#d1c1e9", "TC" = "#241144", "Overlap" = "#6f5990"),
    name = "Category")

print(g)

pdf(file = paste(PATH_results, "PG_counts_by_Condition.pdf", sep=""), width = 4, height = 4)
print(g)
dev.off()

#----------------------------

mat <- df

# reorder colData to match matrix
colData <- colData[colData$sampleID %in% colnames(mat), ]
index   <- match(colnames(mat), colData$sampleID)
colData <- colData[index, ]

head(colnames(mat))
head(colData$sampleID)

Conditions <- unique(colData$Condition)

rows_to_keep <- apply(mat, 1, function(row) {
  any(sapply(Conditions, function(Condition) {
    
    # Subset mat for samples corresponding to the current Condition
    Condition_samples <- colData$sampleID[colData$Condition == Condition & colData$Turbo == "T"]
    Condition_values <- row[match(Condition_samples, colnames(mat))]
    
    # Check if 75% of samples in this Condition have non-NA values
    mean(!is.na(Condition_values)) >= 0.75
  }))
})

# Filter mat to include only rows meeting the condition
mat <- mat[rows_to_keep, ]

# Verify dimensions of the filtered matrix
dim(mat)

# colData MUST match order of mat
head(colnames(mat))
head(colData$sampleID)

colData$ConditionTurbo <- paste0(colData$Condition, "_", colData$Turbo)

design <- model.matrix(~ 0 + ConditionTurbo, data = colData) #select design
head(design)

fit <- lmFit(mat, design)

colnames(fit) #check possible comparisons

#----------------

# Oxaliplatin
contrast_matrix <- makeContrasts( ConditionTurboOx_T -  ConditionTurboOx_TC, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust="BH", number=Inf)

head(results)

sig_de <- results[results$adj.P.Val < 0.05 & results$logFC > 1, , drop = FALSE] #extract significant proteins if there
head(sig_de)

sig_de <- na.omit(sig_de)
dim(sig_de)

# volcano plotting
res       <- as.data.frame(results) 
mutateddf <- mutate(res, Sig=ifelse(res$adj.P.Val<0.05 & abs(res$logFC)>1, "FDR < 0.05 & LFC > 1", ifelse("NS")))
input     <- cbind(gene=rownames(res), mutateddf) 

volc = ggplot(input, aes(logFC, -log10(P.Value))) + geom_point(aes(col=Sig)) +
  #scale_color_manual(values = c("#0D0887FF","#9512A1FF", "grey")) + 
  scale_color_manual(values = c("#9357f1", "grey")) + 
  ggrepel::geom_text_repel(data=subset(input, input$gene %in% rownames(sig_de)),
                           aes(label=gene), size=4, segment.alpha= 0.2, force =2, max.overlaps=16) 
volc <- volc + theme_bw() + theme(aspect.ratio=1)
volc <- volc + theme(legend.position="bottom", axis.text.y = element_text(size= 12, face="bold"), 
                     axis.title.y = element_text(size=14), axis.title.x = element_text(size= 14), 
                     axis.text.x = element_text(size= 12), legend.title=element_text(size=14), 
                     legend.text=element_text(size=14), plot.title=element_text(size=12, hjust = 0.5)) + 
  ggtitle("Oxaliplatin")

print(volc)

pdf(file = paste0(PATH_results, "volcano-ox.pdf"), height = 4, width = 4)
print(volc)
dev.off()

pdf(file = paste0(PATH_results, "volcano-ox_big.pdf"), height = 6, width = 6)
print(volc)
dev.off()

table(sig_de$logFC > 1)
write.csv(results, paste(PATH_results, "DEP-analysis-limma_ox.csv"))
write.csv(sig_de,  paste(PATH_results, "DEP-analysis-limma_ox_sig.csv"))

#---------

# Run hypothesis testing, adjust Condition as needed
contrast_matrix <- makeContrasts( ConditionTurboV_T -  ConditionTurboV_TC, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust="BH", number=Inf)

head(results)

sig_de <- results[results$adj.P.Val < 0.05 & results$logFC > 1, , drop = FALSE] #extract significant proteins if there
head(sig_de)

sig_de <- na.omit(sig_de)
dim(sig_de)

# volcano plotting
res       <- as.data.frame(results) 
mutateddf <- mutate(res, Sig=ifelse(res$adj.P.Val<0.05 & abs(res$logFC)>1, "FDR < 0.05 & LFC > 1", ifelse("NS")))
input     <- cbind(gene=rownames(res), mutateddf) 

volc = ggplot(input, aes(logFC, -log10(P.Value))) + geom_point(aes(col=Sig)) +
  #scale_color_manual(values = c("#0D0887FF","#9512A1FF", "grey")) + 
  scale_color_manual(values = c("#9357f1", "grey")) + 
  ggrepel::geom_text_repel(data=subset(input, input$gene %in% rownames(sig_de)),
                           aes(label=gene), size=4, segment.alpha= 0.2, force =2, max.overlaps=16) 
volc <- volc + theme_bw() + theme(aspect.ratio=1)
volc <- volc + theme(legend.position="bottom", axis.text.y = element_text(size= 12, face="bold"), 
                     axis.title.y = element_text(size=14), axis.title.x = element_text(size= 14), 
                     axis.text.x = element_text(size= 12), legend.title=element_text(size=14), 
                     legend.text=element_text(size=14), plot.title=element_text(size=12, hjust = 0.5)) + 
  ggtitle("Vehicle")

print(volc)

pdf(file = paste0(PATH_results, "volcano-vehicle.pdf"), height = 4, width = 4)
print(volc)
dev.off()

pdf(file = paste0(PATH_results, "volcano-vehicle_big.pdf"), height = 6, width = 6)
print(volc)
dev.off()

table(sig_de$logFC > 1)
write.csv(results, paste(PATH_results, "DEP-analysis-limma_vehicle.csv"))
write.csv(sig_de,  paste(PATH_results, "DEP-analysis-limma_vehicle_sig.csv"))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

dir.create("./output/explants/tissue-enr")
PATH_results = "./output/explants/tissue-enr/"

enrichments <- read.csv("./output/enrichments_75filt.csv", 
                        check.names = FALSE, header = TRUE, row.names = 1)

# filter for ANY enriched gene in the 4 tissue list, instead of a 75% threshold in explants

mat <- df

# reorder colData to match matrix
colData <- colData[colData$sampleID %in% colnames(mat), ]
index   <- match(colnames(mat), colData$sampleID)
colData <- colData[index, ]

Conditions <- unique(colData$Condition)

rows_to_keep 

mat <- mat[rownames(mat) %in% enrichments$Gene, ]

dim(mat)

head(colnames(mat))
head(colData$sampleID)

colData$ConditionTurbo <- paste0(colData$Condition, "_", colData$Turbo)

design <- model.matrix(~ 0 + ConditionTurbo, data = colData) #select design
head(design)

fit <- lmFit(mat, design)

colnames(fit) #check possible comparisons

#----------------

# Oxaliplatin
contrast_matrix <- makeContrasts( ConditionTurboOx_T -  ConditionTurboOx_TC, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust="BH", number=Inf)

head(results)

sig_de <- results[results$adj.P.Val < 0.05 & results$logFC > 1, , drop = FALSE] #extract significant proteins if there
head(sig_de)

sig_de <- na.omit(sig_de)
dim(sig_de)

# volcano plotting
res       <- as.data.frame(results) 
mutateddf <- mutate(res, Sig=ifelse(res$adj.P.Val<0.05 & abs(res$logFC)>1, "FDR < 0.05 & LFC > 1", ifelse("NS")))
input     <- cbind(gene=rownames(res), mutateddf) 

volc = ggplot(input, aes(logFC, -log10(P.Value))) + geom_point(aes(col=Sig)) +
  #scale_color_manual(values = c("#0D0887FF","#9512A1FF", "grey")) + 
  scale_color_manual(values = c("#9357f1", "grey")) + 
  ggrepel::geom_text_repel(data=subset(input, input$gene %in% rownames(sig_de)),
                           aes(label=gene), size=4, segment.alpha= 0.2, force =2, max.overlaps=16) 
volc <- volc + theme_bw() + theme(aspect.ratio=1)
volc <- volc + theme(legend.position="bottom", axis.text.y = element_text(size= 12, face="bold"), 
                     axis.title.y = element_text(size=14), axis.title.x = element_text(size= 14), 
                     axis.text.x = element_text(size= 12), legend.title=element_text(size=14), 
                     legend.text=element_text(size=14), plot.title=element_text(size=12, hjust = 0.5)) + 
  ggtitle("Oxaliplatin")

print(volc)

pdf(file = paste0(PATH_results, "volcano-ox.pdf"), height = 4, width = 4)
print(volc)
dev.off()

pdf(file = paste0(PATH_results, "volcano-ox_big.pdf"), height = 6, width = 6)
print(volc)
dev.off()

table(sig_de$logFC > 1)
write.csv(results, paste(PATH_results, "DEP-analysis-limma_ox.csv"))
write.csv(sig_de,  paste(PATH_results, "DEP-analysis-limma_ox_sig.csv"))

#---------

# Run hypothesis testing, adjust Condition as needed
contrast_matrix <- makeContrasts( ConditionTurboV_T -  ConditionTurboV_TC, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust="BH", number=Inf)

head(results)

sig_de <- results[results$adj.P.Val < 0.05 & results$logFC > 1, , drop = FALSE] #extract significant proteins if there
head(sig_de)

sig_de <- na.omit(sig_de)
dim(sig_de)

# volcano plotting
res       <- as.data.frame(results) 
mutateddf <- mutate(res, Sig=ifelse(res$adj.P.Val<0.05 & abs(res$logFC)>1, "FDR < 0.05 & LFC > 1", ifelse("NS")))
input     <- cbind(gene=rownames(res), mutateddf) 

volc = ggplot(input, aes(logFC, -log10(P.Value))) + geom_point(aes(col=Sig)) +
  #scale_color_manual(values = c("#0D0887FF","#9512A1FF", "grey")) + 
  scale_color_manual(values = c("#9357f1", "grey")) + 
  ggrepel::geom_text_repel(data=subset(input, input$gene %in% rownames(sig_de)),
                           aes(label=gene), size=4, segment.alpha= 0.2, force =2, max.overlaps=16) 
volc <- volc + theme_bw() + theme(aspect.ratio=1)
volc <- volc + theme(legend.position="bottom", axis.text.y = element_text(size= 12, face="bold"), 
                     axis.title.y = element_text(size=14), axis.title.x = element_text(size= 14), 
                     axis.text.x = element_text(size= 12), legend.title=element_text(size=14), 
                     legend.text=element_text(size=14), plot.title=element_text(size=12, hjust = 0.5)) + 
  ggtitle("Vehicle")

print(volc)

pdf(file = paste0(PATH_results, "volcano-vehicle.pdf"), height = 4, width = 4)
print(volc)
dev.off()

pdf(file = paste0(PATH_results, "volcano-vehicle_big.pdf"), height = 6, width = 6)
print(volc)
dev.off()

table(sig_de$logFC > 1)
write.csv(results, paste(PATH_results, "DEP-analysis-limma_vehicle.csv"))
write.csv(sig_de,  paste(PATH_results, "DEP-analysis-limma_vehicle_sig.csv"))
