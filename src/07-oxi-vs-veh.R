library(dplyr)
library(tidyr)
library(limma)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(stringr)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)

PATH_results = "./output/explants/"

df      <- read.csv("./data/explants/combined_genes.csv", row.names = 1, check.names = FALSE)
colData <- read.csv("./data/explants/colData.csv", row.names = 1)

enrichments <- read.csv("./output/explants/tissue-enr/enrichments.csv", row.names = 1)

#-------------------------------------------------------------------------------

head(enrichments)

# Extract neuronal enriched genes (see 06-explant-enrichment.R)
mat <- df[rownames(df) %in% enrichments$Gene, ]

colData <- colData[colData$sampleID %in% colnames(mat), ]
index   <- match(colnames(mat), colData$sampleID)
colData <- colData[index, ]

head(colnames(mat))
head(colData$sampleID)

TC <- colData$sampleID[colData$Turbo == "T"]

mat <- mat[, colnames(mat) %in% TC]
colData <- colData[colData$sampleID %in% TC, ]

design <- model.matrix(~ 0 + Condition, data = colData) #select design
head(design)

fit <- lmFit(mat, design)

colnames(fit) #check possible comparisons

#----------------

# Oxaliplatin
contrast_matrix <- makeContrasts(ConditionOx -  ConditionV, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust="BH", number=Inf)

head(results)

sig_de <- results[results$adj.P.Val < 0.05, , drop = FALSE] #extract significant proteins if there
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
  ggtitle("Oxaliplatin vs Vehicle")

print(volc)

pdf(file = paste0(PATH_results, "volcano.pdf"), height = 4, width = 4)
print(volc)
dev.off()

pdf(file = paste0(PATH_results, "volcano_big.pdf"), height = 6, width = 6)
print(volc)
dev.off()

table(sig_de$adj.P.Val < 0.05)

write.csv(results, paste(PATH_results, "DEP-analysis-limma_drugeffect.csv"))
write.csv(sig_de,  paste(PATH_results, "DEP-analysis-limma_drugeffect_sig.csv"))

#-------------------------------------------------------------------------------

# Biological pathway differences
# NOTE: empty for UP

test_df    <- rownames(sig_de)[sig_de$logFC < -0.5] #adjust for up, down regulated genes
background <- rownames(results) # only test genes that could be tested with limma, not all enriched

ego <- enrichGO(gene          = test_df,
                universe      = background,
                OrgDb         = org.Mm.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

# remove redundancy in the GO terms
ego2 <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min, measure = 'Wang')

pdf(file = paste(PATH_results, "BP_UP-network.pdf", sep=""), width = 8, height = 8)
goplot(ego2)
dev.off()

pdf(file = paste(PATH_results, "BP_UP-barplot.pdf", sep=""), width = 6, height = 3)
mutate(ego2, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")
dev.off()

pdf(file = paste(PATH_results, "BP_UP-upset.pdf", sep=""), width = 15, height = 4)
upsetplot(ego2)
dev.off()

write.csv(ego2, paste(PATH_results, "BP_UP.csv"))


#--------

# Cellular Compartment
# NOTE: empty for UP
ego <- enrichGO(gene          = test_df,
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

pdf(file = paste(PATH_results, "Compartment_DOWN-network.pdf", sep=""), width = 8, height = 8)
goplot(ego2)
dev.off()

pdf(file = paste(PATH_results, "Compartment_DOWN-barplot.pdf", sep=""), width = 6, height = 3)
mutate(ego2, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")
dev.off()

pdf(file = paste(PATH_results, "Compartment_DOWN-upset.pdf", sep=""), width = 15, height = 4)
upsetplot(ego2)
dev.off()

write.csv(ego2, paste(PATH_results, "Compartment_DOWN.csv"))

#--------

# Molecular Function
ego <- enrichGO(gene          = test_df,
                universe      = background,
                OrgDb         = org.Mm.eg.db,
                keyType       = "SYMBOL",
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

# remove redundancy in the GO terms
ego2    <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min, measure = 'Wang')

pdf(file = paste(PATH_results, "MF_UP-network.pdf", sep=""), width = 8, height = 8)
goplot(ego2)
dev.off()

pdf(file = paste(PATH_results, "MF_UP-barplot.pdf", sep=""), width = 6, height = 3)
mutate(ego2, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")
dev.off()

pdf(file = paste(PATH_results, "MF_UP-upset.pdf", sep=""), width = 15, height = 4)
upsetplot(ego2)
dev.off()

write.csv(ego2, paste(PATH_results, "MF_UP.csv"))

#-------------------------------------------------------------------------------
