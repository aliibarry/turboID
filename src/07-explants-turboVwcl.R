results.turbo <- read.csv("./output/explants/ DEP-analysis-limma_drugeffect.csv", row.names = 1)
results.wcl   <- read.csv("./output/explants/wcl/ DEP-analysis-limma_ox.csv", row.names = 1)


turbo <- results.turbo[results.turbo$adj.P.Val < 0.05 & abs(results.turbo$logFC) > 0.5, ]
turbo <- na.omit(turbo)

wcl   <- results.wcl[results.wcl$adj.P.Val < 0.05 & abs(results.wcl$logFC) > 0.5, ]
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
