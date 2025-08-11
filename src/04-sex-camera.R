library(dplyr)
library(tidyr)
library(limma)
library(ggplot2)
library("clusterProfiler")
library(msigdbr)

PATH_results = "./output/sex/"

df          <- read.csv("./data/matrix-for-limma.csv", header = TRUE)
colData     <- read.csv("./data/colData-for-limma.csv", header = TRUE)
enrichments <- read.csv("./output/enrichments.csv", 
                        check.names = FALSE, header = TRUE, row.names = 1)

hallmark_sets <- msigdbr(species = "Mus musculus", category = "H")
hallmark_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)

GO_gene_sets <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
GO_BP_list   <- split(GO_gene_sets$gene_symbol, GO_gene_sets$gs_name)

#-------------------------------------------------------------------------------

head(enrichments)

enrichment <- enrichments[enrichments$Tissue == "DRG",]
df_sub <- df[df$genes %in% enrichment$Gene, ]

head(df_sub)

colData <-  colData[colData$Turbo  == "T", ]
colData <-  colData[colData$Tissue == "DRG", ]

colData$Tissue.Sex <- as.factor(paste(colData$Tissue,colData$Sex, sep = "."))

mat <- df_sub %>% distinct(genes, .keep_all = TRUE) #6109 distinct genes
rownames(mat) <- mat$genes

mat$proteins <- NULL
mat$genes    <- NULL

mat <- mat[,colnames(mat) %in% colData$sampleID]

index   <- match(colnames(mat), colData$sampleID)
colData <- colData[index, ]

head(colData)
mat[1:5, 1:5]

design <- model.matrix(~ 0 + Tissue.Sex, data = colData) #select design

mat <- na.omit(mat)

# mat[is.na(mat)] <- 0

camera_results_HM <- camera(mat, hallmark_list, design)
camera_results_BP <- camera(mat, GO_BP_list , design) 

write.csv(camera_results_BP, paste0(PATH_results, "camera_DRG.csv"), row.names = TRUE)

#-------------------------------------------------------------------------------

colData     <- read.csv("./data/colData-for-limma.csv", header = TRUE)

head(enrichments)

enrichment <- enrichments[enrichments$Tissue == "paw",]
df_sub <- df[df$genes %in% enrichment$Gene, ]

head(df_sub)

colData <-  colData[colData$Turbo  == "T", ]
colData <-  colData[colData$Tissue == "paw", ]

colData$Tissue.Sex <- as.factor(paste(colData$Tissue,colData$Sex, sep = "."))

mat <- df_sub %>% distinct(genes, .keep_all = TRUE) #6109 distinct genes
rownames(mat) <- mat$genes

mat$proteins <- NULL
mat$genes    <- NULL

mat <- mat[,colnames(mat) %in% colData$sampleID]

index   <- match(colnames(mat), colData$sampleID)
colData <- colData[index, ]

head(colData)
mat[1:5, 1:5]

design <- model.matrix(~ 0 + Tissue.Sex, data = colData) #select design

mat[is.na(mat)] <- 0

camera_results_HM <- camera(mat, hallmark_list, design)
camera_results_BP <- camera(mat, GO_BP_list , design) 

write.csv(camera_results_BP, paste0(PATH_results, "camera_paw.csv"), row.names = TRUE)
