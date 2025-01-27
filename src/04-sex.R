library(dplyr)
library(tidyr)
library(limma)
library(ggplot2)
library("clusterProfiler")
library(msigdbr)


PATH_results = "./output/sex/"

df          <- read.csv("./data/matrix-for-limma.csv", header = TRUE)
colData     <- read.csv("./data/colData-for-limma.csv", header = TRUE)
enrichments <- read.csv("./output/enrichments_75filt.csv", 
                        check.names = FALSE, header = TRUE, row.names = 1)

# Combined Turbo + WCL across tissue (ie., all possible proteins that could have been enriched)
background <- readxl::read_excel("./data/JRS_curated/20240826_corBackground_for_GO_WCL-TurboALL.xlsx", col_names = FALSE)
names(background)[names(background) == "...1"] <- "Gene"

#-------------------------------------------------------------------------------

head(enrichments)

enrichments <- enrichments[enrichments$Tissue == "DRG"]
df <- df[df$genes %in% enrichments$Gene, ]

head(df)

colData <-  colData[colData$Turbo  == "T", ]
#colData <-  colData[colData$Tissue == "DRG", ]

colData$Tissue.Sex <- as.factor(paste(colData$Tissue,colData$Sex, sep = "."))

mat <- df %>% distinct(genes, .keep_all = TRUE) #6109 distinct genes
rownames(mat) <- mat$genes

mat$proteins <- NULL
mat$genes    <- NULL

mat <- mat[,colnames(mat) %in% colData$sampleID]

index   <- match(colnames(mat), colData$sampleID)
colData <- colData[index, ]

head(colData)
mat[1:5, 1:5]

design <- model.matrix(~ 0 + Tissue.Sex, data = colData) #select design
head(design)

fit <- lmFit(mat, design)

colnames(fit) #check possible comparisons

# Run hypothesis testing, adjust tissue as needed
contrast_matrix <- makeContrasts( Tissue.SexDRG.m -  Tissue.SexDRG.f, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust="BH", number=Inf)

head(results)

sig_de <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, , drop = FALSE] #extract significant proteins if there
head(sig_de)

sig_de <- na.omit(sig_de)
dim(sig_de)

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
  ggtitle("DRG.sex")

print(volc)

pdf(file = paste0(PATH_results, "volcano-DRG.sex.pdf"), height = 4, width = 4)
print(volc)
dev.off()

pdf(file = paste0(PATH_results, "volcano-DRG.sex_big.pdf"), height = 6, width = 6)
print(volc)
dev.off()

table(sig_de$logFC > 1)
write.csv(results, paste(PATH_results, "DEP-analysis-limma-DRG.sex.csv"))
write.csv(sig_de,  paste(PATH_results, "DEP-analysis-limma-DRG.sex_sig.csv"))

#-------------------------------------------------------------------------------

input <- read.csv("./output/sex/ DEP-analysis-limma-DRG.sex.csv", header = TRUE)
input <- input[c("X", "adj.P.Val", "logFC")]

hallmark_sets = msigdbr(species = "Mus musculus", category = "H")
GO_gene_sets  = msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")

filtered_dge <- input %>%
  dplyr::arrange(dplyr::desc(abs(logFC))) %>%
  dplyr::distinct(X, .keep_all = TRUE)

lfc_vector <- filtered_dge$logFC
names(lfc_vector) <- filtered_dge$X

# Sort log2 fold change values in descending order
lfc_vector <- sort(lfc_vector, decreasing = TRUE)

set.seed(52)

gsea_results <- GSEA(
  geneList = lfc_vector,
  minGSSize = 25, 
  #maxGSSize = 500,
  pvalueCutoff = 0.05, 
  eps = 0, 
  seed = TRUE, 
  pAdjustMethod = "BH",
  TERM2GENE = dplyr::select(hallmark_sets, gs_name, gene_symbol)
)

gsea_result_df <- data.frame(gsea_results@result)
write.csv(gsea_result_df, file = paste0(PATH_results, "GSEA_results_DRG.sex.csv"))

g <- ggplot(gsea_result_df, aes(x=(Description), y=NES, colour=p.adjust, size=setSize))
g <- g + geom_point() + theme_bw() + ggtitle("DRG, TurboID, sex") +
  theme(axis.text.y = element_text(size= 12, colour= "black", hjust = 1), 
        axis.text.x = element_text(size=10, angle = 45, hjust= 1), 
        legend.text = element_text(size=10), 
        axis.title.x = element_blank(),
        plot.title=element_text(size=rel(1), hjust = 1)) +
  theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm")) +  
  labs(y="Enrichment Score", colour="p value", size="Count")

pdf(paste0(PATH_results, "GSEA_dots_hallmark.pdf"), height = 5, width = 4)
print(g)
dev.off()

#----------

gsea_results <- GSEA(
  geneList = lfc_vector,
  minGSSize = 25, 
  #maxGSSize = 500,
  pvalueCutoff = 0.05, 
  eps = 0, 
  seed = TRUE, 
  pAdjustMethod = "BH",
  TERM2GENE = dplyr::select(GO_gene_sets, gs_name, gene_symbol)
)

gsea_result_df <- data.frame(gsea_results@result)
write.csv(gsea_result_df, file = paste0(PATH_results, "GSEA_results_DRG-BP.sex.csv"))

g <- ggplot(gsea_result_df, aes(x=(Description), y=NES, colour=p.adjust, size=setSize))
g <- g + geom_point() + theme_bw() + ggtitle("DRG, TurboID, sex") +
  theme(axis.text.y = element_text(size= 12, colour= "black", hjust = 1), 
        axis.text.x = element_text(size=10, angle = 45, hjust= 1), 
        legend.text = element_text(size=10), 
        axis.title.x = element_blank(),
        plot.title=element_text(size=rel(1), hjust = 1)) +
  theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm")) +  
  labs(y="Enrichment Score", colour="p.adj", size="Count")

pdf(paste0(PATH_results, "GSEA_dots_BP.pdf"), height = 6, width = 4)
print(g)
dev.off()

