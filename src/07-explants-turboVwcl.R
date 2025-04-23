results.turbo <- read.csv("./output/explants/ DEP-analysis-limma_drugeffect.csv", row.names = 1)
results.wcl   <- read.csv("./output/explants/wcl/ DEP-analysis-limma_ox.csv", row.names = 1)


turbo <- results.turbo[results.turbo$adj.P.Val < 0.05 & abs(results.turbo$logFC) > 0, ]
turbo <- na.omit(turbo)

wcl   <- results.wcl[results.wcl$adj.P.Val < 0.05 & abs(results.wcl$logFC) > 0, ]
wcl   <- na.omit(wcl)

merged_data <- merge(turbo, wcl, by = "row.names", all = TRUE, suffixes = c(".turbo",".wcl"))
head(merged_data)

dep_common <- intersect(rownames(turbo), rownames(wcl))
dep_turbo  <- setdiff(rownames(turbo), rownames(wcl))

write.csv(dep_turbo, "./output/dep_turbo-unique.csv")

#-------------------------------------------------------------------------------

results.turbo <- read.csv("./output/explants/ DEP-analysis-limma_drugeffect.csv", row.names = 1)
results.wcl   <- read.csv("./output/explants/wcl/ DEP-analysis-limma_ox.csv", row.names = 1)

turbo <- results.turbo
turbo <- na.omit(turbo)

wcl   <- results.wcl
wcl   <- na.omit(wcl)

merged_data <- merge(turbo, wcl, by = "row.names", all = TRUE, suffixes = c(".turbo",".wcl"))
head(merged_data)

dep_common <- intersect(rownames(turbo), rownames(wcl))
dep_turbo  <- setdiff(rownames(turbo), rownames(wcl))

#-------------------------------------------------------------------------------

library(RRHO2)
library(org.Mm.eg.db)
library(clusterProfiler)

dir.create("./output/explants/RRHO/")
PATH_results = "./output/explants/RRHO/"

# add direction of effect to p value
results.turbo$dde <- -log10(results.turbo$P.Value) * sign(results.turbo$logFC)
results.wcl$dde   <- -log10(results.wcl$P.Value) * sign(results.wcl$logFC)

input.turbo <- data.frame(Genes = rownames(results.turbo),
                          DDE   = results.turbo$dde,
                          stringsAsFactors = FALSE)

input.wcl   <- data.frame(Genes = rownames(results.wcl),
                          DDE   = results.wcl$dde,
                          stringsAsFactors = FALSE)

input.turbo <- na.omit(input.turbo)
input.wcl   <- na.omit(input.wcl)

# filter for overlapping proteins
shared_genes <- intersect(input.turbo$Genes, input.wcl$Genes)
input.turbo  <- input.turbo[input.turbo$Genes %in% shared_genes, ]
input.wcl    <- input.wcl[input.wcl$Genes %in% shared_genes, ]

RRHO_obj <-  RRHO2_initialize(input.turbo, input.wcl, 
                              labels = c("Neuronal", "Bulk"), 
                              log10.ind=TRUE,
                              boundary = 0.025)

RRHO2_heatmap(RRHO_obj)
RRHO2_vennDiagram(RRHO_obj, type="uu")
RRHO2_vennDiagram(RRHO_obj, type="dd")
RRHO2_vennDiagram(RRHO_obj, type="ud")
RRHO2_vennDiagram(RRHO_obj, type="du")

# save figures
pdf(file = paste0(PATH_results, "heatmap.pdf"), height = 6, width = 6)
RRHO2_heatmap(RRHO_obj)
dev.off()

pdf(file = paste0(PATH_results, "concordant_up.pdf"), height = 3, width = 3)
RRHO2_vennDiagram(RRHO_obj, type = "uu")
dev.off()

pdf(file = paste0(PATH_results, "concordant_down.pdf"), height = 3, width = 3)
RRHO2_vennDiagram(RRHO_obj, type = "dd")
dev.off()

pdf(file = paste0(PATH_results, "disconcordant_turbo-up.pdf"), height = 3, width = 3)
RRHO2_vennDiagram(RRHO_obj, type = "ud")
dev.off()

pdf(file = paste0(PATH_results, "disconcordant_turbo-down.pdf"), height = 3, width = 3)
RRHO2_vennDiagram(RRHO_obj, type = "du")
dev.off()

#----------------

# check discordant proteins
du.list    <- RRHO_obj$genelist_du$gene_list_overlap_du
ud.list    <- RRHO_obj$genelist_ud$gene_list_overlap_ud
background <- shared_genes

ego <- enrichGO(gene          = append(du.list,ud.list),
                universe      = background,
                OrgDb         = org.Mm.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

#--

ego.ud <- enrichGO(gene          = ud.list,
                   universe      = background,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = "SYMBOL",
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)
head(ego.ud)
