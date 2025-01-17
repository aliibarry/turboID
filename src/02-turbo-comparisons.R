library(dplyr)
library(tidyr)
library(ggplot2)
library(VennDiagram)

# Compare Type B analyses with original t-test approach

genes <- read.csv("./output/enrichments_75filt.csv", row.names = 1) 
jrs   <- read.csv("./data/JRS_curated/JRS_enriched.csv", header = TRUE)

PATH_results = "./output/enrichments/"

#-------------------------------------------------------------------------------

jrs <- jrs %>%
  separate(Gene, into = paste0("Gene"), sep = ";", remove = TRUE) %>%
  distinct(Gene, Tissue, .keep_all = TRUE)

jrs$analysis <- "t-test"
genes$analysis <- "limma"

merged_genes <- bind_rows(genes, jrs)

tail(genes)

merged_df1 <- merged_genes[merged_genes$Type  == "B", ]
merged_df2 <- merged_genes[merged_genes$Type  == "DEP", ]

merged_df <- bind_rows(merged_df1, merged_df2)

# Box plotting
plot_data <- merged_df %>%
  dplyr::group_by(Tissue, Type, analysis) %>%
  dplyr::summarise(Count = n()) %>%
  dplyr::ungroup()

g <- ggplot(plot_data, aes(x = Tissue, y = Count, fill = analysis)) 
g <- g + geom_bar(stat = "identity", position = "dodge") 
# g <- g + geom_text(aes(label = Count),
#                    position = position_dodge(width = 0.9),  # Align text with bars
#                    vjust = -0.3,                          # Adjust vertical position
#                    color = "black") 
g <- g + labs(title = "Turbo Enrichments, Type B", x = "Tissue", y = "Gene Count") 
g <- g + theme_bw() +
  scale_fill_manual(
    values = c("limma" = "#dfbee1", "t-test" = "#95459b"),
    name = "Analysis Method")

print(g)

pdf(file = paste0(PATH_results, "enrichments_compared.pdf"), width = 5, height = 4)
print(g)
dev.off() 

#----------

tissue_groups <- split(merged_df, merged_df$Tissue)

# Generate Venn diagrams for each Tissue
venn_plots <- list()
for (tissue in names(tissue_groups)) {
  tissue_data <- tissue_groups[[tissue]]
  
  # Create a list of gene sets for each analysis
  venn_list <- split(tissue_data$Gene, tissue_data$analysis)
  
  venn <- venn.diagram(
    x = venn_list,
    category.names = names(venn_list),
    fill = c("#dfbee1", "#95459b"),
    alpha = 0.5,
    cex = 1.0,      
    cat.cex = 0.8,
    margin = 0.1, 
    #main = paste(tissue),
    main.cex = 1.2, 
    filename = NULL   
  )
  
  venn_plots[[tissue]] <- venn
}

grid.newpage()
grid.arrange(
  grobs = venn_plots,
  ncol = 2,  # 2 columns
  nrow = 2   # 2 rows
)

pdf(file = paste0(PATH_results, "enrichments_comparedvenn.pdf"), width = 6, height = 6)
grid.arrange(
  grobs = venn_plots,
  ncol = 2, nrow = 2)
dev.off() 

#-------------------------------------------------------------------------------

merged_genes <- bind_rows(genes, jrs)

tail(genes)

merged_df1 <- merged_genes[merged_genes$Type  == "B", ]
merged_df2 <- merged_genes[merged_genes$Type  == "DEP", ]

merged_df <- bind_rows(merged_df1, merged_df2)

# Box plotting
plot_data <- merged_df %>%
  dplyr::group_by(Tissue, Type, analysis) %>%
  dplyr::summarise(Count = n()) %>%
  dplyr::ungroup()

g <- ggplot(plot_data, aes(x = Tissue, y = Count, fill = analysis)) 
g <- g + geom_bar(stat = "identity", position = "dodge") 
# g <- g + geom_text(aes(label = Count),
#                    position = position_dodge(width = 0.9),  # Align text with bars
#                    vjust = -0.3,                          # Adjust vertical position
#                    color = "black") 
g <- g + labs(title = "Turbo Enrichments, Type B", x = "Tissue", y = "Gene Count") 
g <- g + theme_bw() +
  scale_fill_manual(
    values = c("limma" = "#dfbee1", "t-test" = "#95459b"),
    name = "Analysis Method")

print(g)

#-------------------------------------------------------------------------------

merged_genes <- bind_rows(genes, jrs)

plot_data <- merged_genes %>%
  dplyr::group_by(Tissue, Type, analysis) %>%
  dplyr::summarise(Count = n()) %>%
  dplyr::ungroup()

g <- ggplot(plot_data, aes(x = interaction(analysis, Tissue), y = Count, fill = analysis)) 
g <- g + geom_bar(stat = "identity", position = "stack") 
#g <- g + geom_bar(stat = "identity", aes(group = Type), position = "stack") 
# g <- g + geom_text(aes(label = Count),
#                    position = position_stack(vjust = 0.5),
#                    color = "black") 
g <- g + labs(title = "Turbo Enrichments, A&B", x = " ", y = "Gene Count") 
g <- g + theme_bw() +
  scale_fill_manual(
    values = c("limma" = "#dfbee1", "t-test" = "#95459b"),
    name = "Type")

print(g)

pdf(file = paste0(PATH_results, "enrichments_compared-A-B.pdf"), width = 5, height = 4)
print(g)
dev.off() 
