library(dplyr)
library(tidyr)
library(limma)
library(readxl)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)

PATH_results = "./output/GO/"

df          <- read.csv("./data/matrix-for-limma.csv", header = TRUE)
colData     <- read.csv("./data/colData-for-limma.csv", header = TRUE)
enrichments <- read.csv("./output/enrichments_75filt.csv", 
                        check.names = FALSE, header = TRUE, row.names = 1)

# Combined Turbo + WCL across tissue (ie., all possible proteins that could have been enriched)
background <- read_excel("./data/JRS_curated/20240826_corBackground_for_GO_WCL-TurboALL.xlsx", col_names = FALSE)
names(background)[names(background) == "...1"] <- "Gene"

#-------------------------------------------------------------------------------

test_df <- enrichments$Gene[enrichments$Tissue == "DRG"] 

ego <- enrichGO(gene          = test_df,
                universe      = background$Gene,
                OrgDb         = org.Mm.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
#head(ego)
goplot(ego)
enrichMap(ego)

pdf(file = paste(PATH_results, "GO_paw-network.pdf", sep=""), width = 8, height = 8)
goplot(ego2)
dev.off()

pdf(file = paste(PATH_results, "GO_paw-barplot.pdf", sep=""), width = 6, height = 3)
mutate(ego, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")
dev.off()

pdf(file = paste(PATH_results, "GO-paw-upset.pdf", sep=""), width = 15, height = 4)
upsetplot(ego)
dev.off()

write.csv(ego, paste(PATH_results, "GO_paw.csv"))

#-------------------------------------------------------------------------------

library(GOSemSim)

# dir.create("./output/GOSemSim")
PATH_results = "./output/GO/GOSemSim/"

test_df <- enrichments$Gene[enrichments$Tissue == "LSC"] 

ego <- enrichGO(gene          = test_df,
                universe      = background$Gene,
                OrgDb         = org.Mm.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01,
                readable      = TRUE)

# remove redundancy in the GO terms
ego2    <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min, measure = 'Wang')

results <- ego2@result[ego2@result$p.adjust < 0.05, ]
results$GeneRatio <- sapply(strsplit(results$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
#results <- results[results$GeneRatio > 0.025, ] #terms with 2.5% of enriched list

dim(results)

#---------

pdf(file = paste(PATH_results, "GO_slim_LSC-network.pdf", sep=""), width = 10, height = 8)
goplot(ego2, showCategory = 10)
dev.off()

results$Description <- factor(results$Description, levels = results$Description[order(results$GeneRatio)])

# Plot top (by ratio overlap to enriched list)
g <- ggplot(results[1:25,], aes(x = Description, y = 0, colour = -log10(p.adjust), size = Count))
g <- g + geom_point() 
g <- g + ggtitle("") + theme_minimal() +  
  theme(
    axis.title = element_blank(),      
    axis.text.y = element_blank(),       
    axis.ticks.y = element_blank(),       
    axis.line = element_blank(),         
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5), 
    axis.text.x = element_text(angle = 45, hjust = 1),  
    plot.margin = margin(0, 0, 0, 0)) 
g <- g + scale_y_continuous(expand = c(0, 0), limits = c(0, 0))  # Keep dots strictly at y = 0

print(g)

pdf(paste0(PATH_results, "GO_slim_LSC.pdf"), height = 3, width = 7)
print(g)
dev.off()

write.csv(ego2, paste(PATH_results, "GO_slim_LSC.csv"))

