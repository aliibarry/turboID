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

pdf(file = paste(PATH_results, "GO_DRG-network.pdf", sep=""), width = 8, height = 8)
goplot(ego)
dev.off()

pdf(file = paste(PATH_results, "GO_DRG-barplot.pdf", sep=""), width = 6, height = 3)
mutate(ego, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")
dev.off()

pdf(file = paste(PATH_results, "GO-DRG-upset.pdf", sep=""), width = 15, height = 4)
upsetplot(ego)
dev.off()

write.csv(ego, paste(PATH_results, "GO_DRG.csv"))

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

#-------------------------------------------------------------------------------

# Molecular Function & Cellular Compartments
PATH_results = "./output/GO/"

test_df <- enrichments$Gene[enrichments$Tissue == "paw"] 

ego <- enrichGO(gene          = test_df,
                universe      = background$Gene,
                OrgDb         = org.Mm.eg.db,
                keyType       = "SYMBOL",
                ont           = "CC", #CC or MF
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01,
                readable      = TRUE)

# remove redundancy in the GO terms
ego2    <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min, measure = 'Wang')

pdf(file = paste(PATH_results, "Compartment_paw-network.pdf", sep=""), width = 8, height = 8)
goplot(ego2)
dev.off()

pdf(file = paste(PATH_results, "Compartment_paw-barplot.pdf", sep=""), width = 6, height = 3)
mutate(ego2, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")
dev.off()

pdf(file = paste(PATH_results, "Compartment_paw-upset.pdf", sep=""), width = 15, height = 4)
upsetplot(ego2)
dev.off()

write.csv(ego2, paste(PATH_results, "Compartment_paw.csv"))

#-------------------------------------------------------------------------------

dir.create("./output/GO/unique-enrichments/")
PATH_results = "./output/GO/unique-enrichments/"

# Keep only genes exclusive to tissue of interest
exclusive_genes <- enrichments %>%
  group_by(Gene) %>%    
  filter(all(Tissue == "paw")) %>%           
  pull(Gene)                                   

# Result
ego <- enrichGO(gene          = exclusive_genes,
                universe      = background$Gene,
                OrgDb         = org.Mm.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP", #BP, CC or MF
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01,
                readable      = TRUE)

ego2    <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min, measure = 'Wang')

write.csv(ego2, paste(PATH_results, "GO_paw_unique.csv"))

pdf(file = paste(PATH_results, "GO_paw-barplot.pdf", sep=""), width = 6, height = 3)
mutate(ego2, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")
dev.off()

pdf(file = paste(PATH_results, "GO_paw-upset.pdf", sep=""), width = 8, height = 4)
upsetplot(ego2)
dev.off()

#-------------------------------------------------------------------------------

PATH_results = "./output/GO/unique-enrichments/"

# Keep only genes exclusive to tissue of interest
exclusive_genes <- enrichments %>%
  group_by(Gene) %>%                     
  filter(n_distinct(Tissue) == n_distinct(enrichments$Tissue)) %>% 
  ungroup() %>%                              
  distinct(Gene)                                                

# Result
ego <- enrichGO(gene          = exclusive_genes$Gene,
                universe      = background$Gene,
                OrgDb         = org.Mm.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP", #BP, CC or MF
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01,
                readable      = TRUE)

ego2    <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min, measure = 'Wang')

write.csv(ego2, paste(PATH_results, "GO_all_unique.csv"))

pdf(file = paste(PATH_results, "GO_all-barplot.pdf", sep=""), width = 6, height = 3)
mutate(ego2, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")
dev.off()

pdf(file = paste(PATH_results, "GO_all-upset.pdf", sep=""), width = 8, height = 4)
upsetplot(ego2)
dev.off()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

PATH_results = "./output/receptors/"

wcl  <- read.csv("./data/JRS_curated/WCL_background.csv", check.names = FALSE)
full <- read_excel("./data/JRS_curated/20240826_corBackground_for_GO_WCL-TurboALL.xlsx", col_names = FALSE)
names(full)[names(full) == "...1"] <- "Gene"
full <- as.data.frame(full)

turbo_enr <- read.csv("./output/enrichments_75filt.csv", check.names = FALSE, header = TRUE, row.names = 1)
turbo_enr <- turbo_enr[, colnames(turbo_enr) %in% c("Gene"), drop = FALSE]

turbo_gen <- read.csv("./data/matrix-for-limma.csv", header = TRUE)
turbo_gen$Gene <- turbo_gen$genes
turbo_gen <- turbo_gen[, colnames(turbo_gen) %in% c("Gene"), drop = FALSE]

wcl$Source       <- "WCL"
turbo_gen$Source <- "Turbo all"
turbo_enr$Source <- "Turbo enriched"

turbo_gen <- as.data.frame(turbo_gen) %>% distinct() 
turbo_enr <- as.data.frame(turbo_enr) %>% distinct() 

df <- rbind(wcl, turbo_enr, turbo_gen)

table(df$Source)

#-------------------------------------------------------------------------------

pg_types <- readxl::read_excel("./data/interactome_list_v3.1_large.xlsx")

rec_types <- data.frame(symbol = pg_types$Rec_symbol, type = pg_types$Rec_type)
lig_types <- data.frame(symbol = pg_types$Lig_symbol, type = pg_types$Lig_type)

pg_types <- rbind(rec_types, lig_types)

head(pg_types)

# Convert interactome list to mouse identifiers using biomart
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

genelist <- pg_types$symbol #extract symbols
genelist <- convertHumanGeneList(genelist)

pg_types <- merge(genelist, pg_types, 
                     by.x = "HGNC.symbol", 
                     by.y = "symbol", 
                     all.x = TRUE)

pg_types <- pg_types %>% distinct() # remove duplicates

head(pg_types)

#-------------------------------------------------------------------------------

type_wcl   <- pg_types[pg_types$MGI.symbol %in% wcl$Gene, ]
type_t.enr <- pg_types[pg_types$MGI.symbol %in% turbo_enr$Gene, ]
type_t.gen <- pg_types[pg_types$MGI.symbol %in% turbo_gen$Gene, ]

type_wcl   <- as.data.frame(table(type_wcl$type))
type_t.enr <- as.data.frame(table(type_t.enr$type))
type_t.gen <- as.data.frame(table(type_t.gen$type))

table(df$Source)
type_wcl   <- type_wcl %>% mutate(Proportion = Freq / 9876)
type_t.enr <- type_t.enr %>% mutate(Proportion = Freq / 4684)
type_t.gen <- type_t.gen %>% mutate(Proportion = Freq / 6109)

type_wcl$Source   <- "WCL"
type_t.gen$Source <- "Turbo all"
type_t.enr$Source <- "Turbo enriched"

combined_counts  <- rbind(type_wcl, type_t.enr, type_t.gen)
head(combined_counts)

g <- ggplot(combined_counts, aes(x = Var1, y = Proportion, fill = Source)) 
g <- g + geom_bar(stat = "identity", position = "dodge") 
g <- g + theme_bw() +
  labs(title = "Receptor prop within each dataset",
       x = "",
       y = "Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(g)

pdf(paste0(PATH_results, "receptor_proportions.pdf"), height = 5, width = 12)
print(g)
dev.off()

#-------------------------------------------------------------------------------

# Look at proteins only seen in the turboID condition
table(df$Source)

PATH_results = "./output/GO/"
background <- read_excel("./data/JRS_curated/20240826_corBackground_for_GO_WCL-TurboALL.xlsx", col_names = FALSE)
names(background)[names(background) == "...1"] <- "Gene"

test_df <- setdiff(unique(df[df$Source == "Turbo all", "Gene"]), unique(df[df$Source == "WCL", "Gene"]))

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

pdf(file = paste(PATH_results, "TurboAll-network.pdf", sep=""), width = 8, height = 8)
goplot(ego2)
dev.off()

pdf(file = paste(PATH_results, "TurboAll-barplot.pdf", sep=""), width = 6, height = 3)
mutate(ego2, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")
dev.off()

pdf(file = paste(PATH_results, "TurboAll-upset.pdf", sep=""), width = 5, height = 4)
upsetplot(ego2)
dev.off()

write.csv(ego2, paste(PATH_results, "TurboAll.csv"))

#-------------------------
