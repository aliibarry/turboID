library(GOSemSim)
library(org.Mm.eg.db)
library(biomaRt)

whole <- read.csv("./output/explants/wcl/ DEP-analysis-limma_ox.csv", row.names = 1)
turbo <- read.csv("./output/explants/ DEP-analysis-limma_drugeffect.csv", row.names = 1)

whole <- na.omit(whole)
turbo <- na.omit(turbo)

d <- GOSemSim::godata('org.Mm.eg.db', ont="CC", computeIC=FALSE)
cluster1 <- rownames(turbo)[turbo$adj.P.Val < 0.05 & abs(turbo$logFC) > 0.5]
cluster2 <- rownames(whole)[whole$adj.P.Val < 0.05 & abs(whole$logFC) > 0.5]

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
cluster1 <- getBM(attributes = c("external_gene_name", "entrezgene_id"),
                    filters = "external_gene_name",
                    values = cluster1,
                    mart = mart)

cluster2 <- getBM(attributes = c("external_gene_name", "entrezgene_id"),
                    filters = "external_gene_name",
                    values = cluster2,
                    mart = mart)

# Calculate similarities between GO pathways for sig_de lists
results_BP <- GOSemSim::clusterSim(cluster1$entrezgene_id, cluster2$entrezgene_id, semData=d, measure="Wang")
write.csv(results_BP, "./output/explants/oxaliplatin/clusterSim_BP.csv")

results_CC <- GOSemSim::clusterSim(cluster1$entrezgene_id, cluster2$entrezgene_id, semData=d, measure="Wang")
write.csv(results_CC, "./output/explants/oxaliplatin/clusterSim_CC.csv")

#-------------------------------------------------------------------------------

library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(ggplot2)
PATH_results = "./output/explants/oxaliplatin/"

#test_df    <- rownames(turbo)[turbo$adj.P.Val < 0.05 & abs(turbo$logFC) > 0.5]
test_df    <- rownames(turbo)[turbo$adj.P.Val < 0.05 & turbo$logFC < -0.5]
background <- rownames(turbo) # only test genes that could be tested with limma, not all enriched

ego_tID1 <- enrichGO(gene     = test_df,
                universe      = background,
                OrgDb         = org.Mm.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

ego_tID2 <- clusterProfiler::simplify(ego_tID1, cutoff=0.7, by="p.adjust", select_fun=min, measure = 'Wang')


fold_tID <- setNames(turbo$logFC, rownames(turbo))
p1 <- cnetplot(ego_tID2, foldChange=fold_changes)

p1 <- cnetplot(ego_tID2, node_label="category", layout = "star", 
               color.params = list(edge = TRUE, category = "black")) 
p2 <- cnetplot(ego_tID2, node_label="gene", layout = "star", 
               color.params = list(foldChange = fold_tID, 
                                   #edge = TRUE, 
                                   category = "black", 
                                   gene ="black")) 

p1 + p2

#-------------------

#test_df    <- rownames(whole)[whole$adj.P.Val < 0.05 & abs(whole$logFC) > 0.5]
test_df    <- rownames(whole)[whole$adj.P.Val < 0.05 & whole$logFC < -0.5]
background <- rownames(whole) # only test genes that could be tested with limma, not all enriched

ego_wcl1 <- enrichGO(gene      = test_df,
                universe      = background,
                OrgDb         = org.Mm.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

ego_wcl2 <- clusterProfiler::simplify(ego_wcl1, cutoff=0.7, by="p.adjust", select_fun=min, measure = 'Wang')

fold_wcl <- setNames(whole$logFC, rownames(whole))

p3 <- cnetplot(ego_wcl2, node_label="category", layout = "star", 
               color.params = list(edge = TRUE, category = "black")) 
p4 <- cnetplot(ego_wcl2, node_label="gene", layout = "star", 
               color.params = list(foldChange = fold_wcl, 
                                   #edge = TRUE, 
                                   category = "black", 
                                   gene ="black")) 
p3 + p4


cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4], rel_widths=c(.8, .8, 1.2))

p5 <- cnetplot(ego_tID2, node_label="category", showCategory = 25,
               color.params = list(edge = TRUE, category = "black")) 

p6 <- cnetplot(ego_wcl2, node_label="category", showCategory = 25,
               color.params = list(edge = TRUE, category = "black")) 

p5 + p6

edo_wcL <- pairwise_termsim(ego_wcl1)
emapplot(edo_wcL)

edo_tid <- pairwise_termsim(ego_tID1)
emapplot(edo_tid)
