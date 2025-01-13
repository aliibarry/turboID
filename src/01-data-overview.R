library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)

dir.create("./output/QC")
PATH_results = "./output/QC/"

# extract all sheets from an excel document
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

#-------------------------------------------------------------------------------

overview   <- read_excel_allsheets("./data/JRS_curated/20240527overview_TurboID_March24.xlsx")
background <- read_excel_allsheets("./data/JRS_curated/20240826_corBackground_for_GO_WCL-TurboALL.xlsx")
paw_rerun  <- read_excel_allsheets("./data/JRS_curated/20241219_rerun_paw_DIANNquant.xlsx")

#-------------------------------------------------------------------------------

df <- read.csv("./data/DIANN_R_output/allTissues_T-TC_Com_gg.csv", check.names = FALSE)

column_names <- colnames(df)[-1]
metadata <- do.call(rbind, strsplit(column_names, "_"))

metadata <- as.data.frame(metadata)
colnames(metadata) <- c("Sample", "Tissue", "Turbo", "Sex", "ID") #clarify TC/T and ID from Julia
metadata$sampleID <- column_names

# Reshape the dataframe to long format
df_long <- df %>%
  pivot_longer(cols = -GeneID, names_to = "sampleID", values_to = "Expression") %>%
  left_join(metadata, by = "sampleID")

summed_expression <- df_long %>%
  group_by(sampleID, Tissue, Turbo) %>%
  summarize(SumCount = sum(Expression, na.rm = TRUE), .groups = "drop") %>%
  mutate(Tissue = factor(Tissue, levels = c("LSC", "DRG", "SCN", "paw")))

mean_intensity <- df_long %>%
  group_by(sampleID, Tissue, Turbo) %>%
  summarize(MeanIntensity = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
  mutate(Tissue = factor(Tissue, levels = c("LSC", "DRG", "SCN", "paw")))

detected_counts <- df_long %>%
  group_by(sampleID, Tissue, Turbo) %>%
  summarize(DetectedCount = sum(!is.na(Expression)), .groups = "drop") %>%
  mutate(Tissue = factor(Tissue, levels = c("LSC", "DRG", "SCN", "paw")))

#-------

g1 <- ggplot(summed_expression, aes(x = Tissue, y = SumCount, fill = Turbo)) 
g1 <- g1 + geom_boxplot() 
g1 <- g1 + theme_classic() +
  labs(title = "", 
       x = "", 
       y = "Summed Expression")

print(g1)

#------

g2 <- ggplot(mean_intensity, aes(x = Tissue, y = MeanIntensity, fill = Turbo)) 
g2 <- g2 + geom_boxplot() 
g2 <- g2 + theme_classic() +
  labs(title = "", 
       x = "", 
       y = "Mean Intensity")

print(g2)

#-------

g3 <- ggplot(detected_counts, aes(x = Tissue, y = DetectedCount, fill = Turbo)) 
g3 <- g3 + geom_boxplot() 
g3 <- g3 + theme_classic() +
  labs(title = "", 
       x = "", 
       y = "PG Count")

print(g3)


#------

pdf(file = paste(PATH_results, "expression_sum.pdf", sep="/"), width = 5, height = 4)
print(g1)
dev.off()

pdf(file = paste(PATH_results, "expression_mean.pdf", sep="/"), width = 5, height = 4)
print(g2)
dev.off()

pdf(file = paste(PATH_results, "PG_counts.pdf", sep="/"), width = 5, height = 4)
print(g3)
dev.off()

#-------------------------------------------------------------------------------



