library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggbiplot)
library(ggrepel)
library(viridis)
library(circlize)
library(ComplexHeatmap)
library(matrixStats)
library(gridExtra)
library(stringr)

dir.create("./output/explants/QC")
PATH_results = "./output/explants/QC/"

df      <- read.csv("./data/explants/combined_genes.csv", row.names = 1, check.names = FALSE)
colData <- read.csv("./data/explants/colData.csv", row.names = 1)

#-------------------------------------------------------------------------------