library(Rsubread)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(genefilter)
library(viridis)
library(readr)
library(gridtext)
require(ComplexHeatmap)
library(circlize)

# setting global variables

c_labels <-
  c("old leaf 1", "old leaf 2", "old leaf 3", "young leaf 1", "young leaf 2", "young leaf 3", "root 1", "root 2", "root 3", "berry 1", "berry 2", "berry 3")
cluster_1 <-
  c(
    "g3303",
    "g3304",
    "g3305",
    "g3306",
    "FUN_112423",
    "g3307",
    "g3308",
    "g3310",
    "g3311",
    "g3312b",
    "g3312",
    "g3313",
    "g3314",
    "g3315",
    "g3318",
    "g3320"
    )
labels_1 <-
  c(
    "CYP87G1-1",
    "24ISO1",
    "SULF2",
    "CYP749B2-1",
    "CYP88C10-1",
    "SDH1",
    "ACT1",
    "CYP88C7-1",
    "CYP88C10-2",
    "CYP87G1-2",
    "CYP87G1-3",
    "24ISO2",
    "SDH2",
    "CYP88C7-2",
    "ACT2",
    "SULF1"
  )
cluster_2 <-
  c(
    "g18181",
    "g18182",
    "g18183",
    "g18184",
    "g18185",
    "g18187",
    "g18188",
    "g18189",
    "g18191",
    "g18192",
    "g18193",
    "g18194",
    "g18195"
  )

labels_2 <-
  c(
    "CYP749B2-5",
    "OGD7",
    "CYP749B2-4",
    "CYP749B2-3",
    "CYP88C8-2",
    "OGD6",
    "OGD5",
    "OGD4",
    "CYP88C7-3",
    "CYP88C8-1",
    "OGD3",
    "OGD2",
    "OGD1"
  )

#Read counts from file
all_counts <- as.matrix(read.delim('normalized_counts_final.txt', header = TRUE, row.names = 1))
head(all_counts)

#Choose only clusters from normalized_counts_final and relabel 
mat <- all_counts[c(cluster_1, cluster_2), ]
row.names(mat) <- c(labels_1, labels_2)
mat_means <- mat - rowMeans(mat)

svg(file = "clusters_functional_centered_updated_og.svg")
#any value over 5 appears red and any value under -5 appears blue
col_fun = colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))
ht_m = Heatmap(
  mat_means,
  name = "log2",
  col = col_fun,
  cluster_rows = F,
  cluster_columns = F,
  row_split = data.frame(rep(c(
    "cluster 1", "cluster 2"
  ), c(16, 13))),
  column_labels = gt_render(c_labels, padding = unit(c(5, 0, 0, 0), "pt")),
)
draw(ht_m)
dev.off()

# CLEAN UP #################################################
# Clear environment
rm(list = ls())

# Clear mind :)
