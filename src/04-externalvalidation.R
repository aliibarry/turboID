library(readxl)
library(biomaRt)
library(VennDiagram)
library(gridExtra)

library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)

# multiple sheets, from Rawat et al. 2025, 10.1097/j.pain.0000000000003743
# paw_ext <- read_excel("./data/published/pain_2025_06_19_sendtner_1_sdc2.xlsx")

paw_ext  <- read.csv("./data/published/rawat2025_mousepawenrichment.csv")
skin_ext <- read.csv("./data/published/rawat2025_humanskinenrichment.csv")

enrichments <- read.csv("./output/enrichments.csv")

PATH_results = "./output/comparisons/"
#-------------------------------------------------------------------------------

# Convert to mouse symbols for comparisons to TurboID
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

# mouse paw samples
genelist <- paw_ext$Genes #extract symbols
genelist <- convertHumanGeneList(genelist)

paw_ext <- merge(genelist, paw_ext, 
                  by.x = "HGNC.symbol", 
                  by.y = "Genes", 
                  all.x = TRUE)

#human skin biopsy
genelist <- skin_ext$Genes #extract symbols
genelist <- convertHumanGeneList(genelist)

skin_ext <- merge(genelist, skin_ext, 
                 by.x = "HGNC.symbol", 
                 by.y = "Genes", 
                 all.x = TRUE)

skin_ext <- skin_ext

#-------------------------------------------------------------------------------

paw_turbo <- enrichments$Gene[enrichments$Tissue == "paw"]

head(paw_ext)
head(paw_turbo)

set1.1 <- unique(paw_ext$MGI.symbol) # anything detected
set1.2 <- unique(paw_ext$MGI.symbol[paw_ext$q.value < 0.05 & paw_ext$log2FC.S2.P1. > 0]) #enriched
set2 <- paw_turbo #turbo enriched in peripheral terminals

turbo_unique_mouse <- setdiff(set2, set1.2)

# Create Venn diagram
venn.plot1 <- venn.diagram(
  x = list(
    paw_ext = set1.1,
    paw_turbo = set2
  ),
  category.names = c("Murine Skin", "Turbo"),
  filename = NULL,   # Keep in R instead of writing to a file
  fill = c("skyblue", "orange"),
  alpha = 0.5,
  cex = 2,
  fontface = "bold",
  cat.cex = 1.5,
  cat.fontface = "bold"
)

venn.plot2 <- venn.diagram(
  x = list(
    paw_ext = set1.2,
    paw_turbo = set2
  ),
  category.names = c("Murine Skin", "Turbo"),
  filename = NULL,   # Keep in R instead of writing to a file
  fill = c("skyblue", "orange"),
  alpha = 0.5,
  cex = 2,
  fontface = "bold",
  cat.cex = 1.5,
  cat.fontface = "bold"
)

# Plot to R graphics window
grid.newpage()
grid.draw(venn.plot1)

pdf(file = paste0(PATH_results, "/rawat2025_mouse_all.pdf"), height = 4, width = 4)
grid.draw(venn.plot1)
dev.off()

pdf(file = paste0(PATH_results, "/rawat2025_mouse_enriched.pdf"), height = 4, width = 4)
grid.draw(venn.plot2)
dev.off()

#-------------------------------------------------------------------------------

set1.1 <- unique(skin_ext$MGI.symbol)
set1.2 <- unique(skin_ext$MGI.symbol[skin_ext$q.value < 0.05 & skin_ext$log2FC.S2.P1. > 0])
set2 <- paw_turbo

turbo_unique_human <- setdiff(set2, set1.2)

# Create Venn diagram
venn.plot1 <- venn.diagram(
  x = list(
    paw_ext = set1.1,
    paw_turbo = set2
  ),
  category.names = c("Human Skin", "Turbo"),
  filename = NULL,   # Keep in R instead of writing to a file
  fill = c("skyblue", "orange"),
  alpha = 0.5,
  cex = 2,
  fontface = "bold",
  cat.cex = 1.5,
  cat.fontface = "bold"
)

venn.plot2 <- venn.diagram(
  x = list(
    paw_ext = set1.2,
    paw_turbo = set2
  ),
  category.names = c("Human Skin", "Turbo"),
  filename = NULL,   # Keep in R instead of writing to a file
  fill = c("skyblue", "orange"),
  alpha = 0.5,
  cex = 2,
  fontface = "bold",
  cat.cex = 1.5,
  cat.fontface = "bold"
)

# Plot to R graphics window
grid.newpage()
grid.draw(venn.plot1)

grid.newpage()
grid.draw(venn.plot2)

pdf(file = paste0(PATH_results, "/rawat2025_human_all.pdf"), height = 4, width = 4)
grid.draw(venn.plot1)
dev.off()

pdf(file = paste0(PATH_results, "/rawat2025_human_enriched.pdf"), height = 4, width = 4)
grid.draw(venn.plot2)
dev.off()

#-------------------------------------------------------------------------------

set1.1 <- unique(paw_ext$MGI.symbol) # anything detected
set1.2 <- unique(paw_ext$MGI.symbol[paw_ext$q.value < 0.05 & paw_ext$log2FC.S2.P1. > 0]) #enriched
set2 <- paw_turbo #turbo enriched in peripheral terminals

turbo_unique_mouse <- setdiff(set2, set1.2)

head(turbo_unique_mouse)
head(turbo_unique_human)

ego <- enrichGO(gene          = turbo_unique_mouse, #intersect(set2, set1.2),
                universe      = append(set1.2, set2), #'set1.2, #'
                OrgDb         = org.Mm.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01,
                readable      = TRUE)

# remove redundancy in the GO terms
ego2    <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min, measure = 'Wang')

write.csv(ego2, "./output/comparisons/TurboVsRawatPaw.non-overlap.BP.csv")

pdf(file = paste(PATH_results, "TurboVsRawatPawintersect_CC-network.pdf", sep=""), width = 8, height = 8)
goplot(ego2)
dev.off()

pdf(file = paste(PATH_results, "TurboVsRawatPawintersect_CC-barplot.pdf", sep=""), width = 6, height = 3)
mutate(ego2, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")
dev.off()

pdf(file = paste(PATH_results, "TurboVsRawatPawintersect_CC-upset.pdf", sep=""), width = 8, height = 4)
upsetplot(ego2)
dev.off()

#-------------------------------------------------------------------------------

set1.1 <- unique(skin_ext$MGI.symbol) # anything detected
set1.2 <- unique(skin_ext$MGI.symbol[skin_ext$q.value < 0.05 & skin_ext$log2FC.S2.P1. > 0]) #enriched
set2 <- paw_turbo #turbo enriched in peripheral terminals

turbo_unique_human <- setdiff(set2, set1.2)

head(turbo_unique_mouse)
head(turbo_unique_human)

ego <- enrichGO(gene          = intersect(set2, set1.2), #set2, #
                universe      = set1.2, # append(set1.2, set2), #
                OrgDb         = org.Mm.eg.db,
                keyType       = "SYMBOL",
                ont           = "CC", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.1,
                qvalueCutoff  = 0.1,
                readable      = TRUE)

# remove redundancy in the GO terms
ego2    <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min, measure = 'Wang')

write.csv(ego2, "./output/comparisons/TurboVsRawatSkin.overlap.CC.csv")

pdf(file = paste(PATH_results, "TurboVsRawatSkinintersect_CC-network.pdf", sep=""), width = 8, height = 8)
goplot(ego2)
dev.off()

pdf(file = paste(PATH_results, "TurboVsRawatSkinintersect_CC-barplot.pdf", sep=""), width = 6, height = 3)
mutate(ego2, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")
dev.off()

pdf(file = paste(PATH_results, "TurboVsRawatSkinintersect_CC-upset.pdf", sep=""), width = 8, height = 4)
upsetplot(ego2)
dev.off()

#-------------------------------------------------------------------------------

# SC synaptosomes (full SC), SOD1 study (Casañas &  Montesinos, 2022)
# https://doi.org/10.1016/j.mcn.2022.103792

library(biomaRt)

sc_ext <- read.csv("./data/published/1-s2.0-S1044743122000987-mmc3.csv")
mart   <- useMart("ensembl", dataset = "mmusculus_gene_ensembl") # for human

mapping <- getBM(
  attributes = c("uniprotswissprot", "mgi_symbol"),
  filters = "uniprotswissprot",
  values = sc_ext$Accession,
  mart = mart
)

mapping

#-----------------

set1 <- mapping$mgi_symbol #SC synaptosome
set2 <- enrichments$Gene[enrichments$Tissue == "LSC"]

turbo_unique_sc <- setdiff(set2, set1)

ego <- enrichGO(gene          = intersect(set2, set1), #set2, #
                universe      = set1, # append(set1.2, set2), #
                OrgDb         = org.Mm.eg.db,
                keyType       = "SYMBOL",
                ont           = "CC", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01,
                readable      = TRUE)

# remove redundancy in the GO terms
ego2    <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min, measure = 'Wang')

write.csv(ego2, "./output/comparisons/TurboVsSC.overlap.CC.csv")

pdf(file = paste(PATH_results, "TurboVsSCintersect_CC-network.pdf", sep=""), width = 8, height = 8)
goplot(ego2)
dev.off()

pdf(file = paste(PATH_results, "TurboVsSCintersect_CC-barplot.pdf", sep=""), width = 6, height = 3)
mutate(ego2, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")
dev.off()

pdf(file = paste(PATH_results, "TurboVsSCintersect_CC-upset.pdf", sep=""), width = 8, height = 4)
upsetplot(ego2)
dev.off()



