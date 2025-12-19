library(Rsubread)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(genefilter)
library(viridis)
library(readr)
library(gridtext)
require(ComplexHeatmap)

# setting global variables

setwd("working_directory_filepath")
c_labels <-
  c("leaf 1", "leaf 2", "leaf 3", "young leaf 1", "young leaf 2", "young leaf 3", "root 1", "root 2", "root 3", "berry 1", "berry 2", "berry 3")
cluster_1 <-
  c(
    "g3303",
    "g3304",
    "g3305",
    "g3306",
    "FUN_112423",
    "g3307",
    "g3308",
    "g3309",
    "g3310",
    "g3311",
    "g3312b",
    "g3312",
    "g3313",
    "g3314",
    "g3315",
    "g3316_g3317",
    "g3318",
    "g3319",
    "g3320"
    )
labels_1 <-
  c(
    "CYP87G1_1",
    "24ISO_1",
    "SULF_2",
    "CYP749B2_1",
    "CYP88C10_1",
    "SDH_1",
    "AT_1",
    "ALP1-like",
    "CYP88C7_1",
    "CYP88C10_2",
    "CYP87G1_2",
    "CYP87G1_3",
    "24ISO_2",
    "SDH_2",
    "CYP88C7_2",
    "psuedo_CYP749B2",
    "AT_2",
    "zinc-finger-like",
    "SULF_1"
  )

cluster_2 <-
  c(
    "g18181",
    "g18182",
    "g18183",
    "g18184",
    "g18185",
    "g18186",
    "g18187",
    "g18188",
    "g18189",
    "g18190",
    "g18191",
    "g18192",
    "g18193",
    "g18194",
    "g18195",
    "g18196",
    "g18197"
  )
labels_2 <-
  c(
    "CYP749B2_5",
    "OGD7",
    "CYP749B2_4",
    "CYP749B2_3",
    "CYP88C8_2",
    "pseudo_CYP749B2",
    "OGD6",
    "OGD5",
    "OGD4",
    "pseudo_OGD",
    "CYP88C7_3",
    "CYP88C8_1",
    "OGD3",
    "OGD2",
    "OGD1",
    "pseudo_24ISO",
    "psuedo_CYP87G1"
  )

# format featureCounts output into usable count matrices

counts_tot <- read.table("featurecounts_multimapping.txt", row.names = 1, header = T)
counts <- counts_tot[,c(6:17)]
counts_r <- round(counts, digits=0) #rounded because I used fractional multimapping mode in featurecounts
colnames(counts_r) <- c("wsleaf_1", "wsleaf_2", "wsleaf_3", "wsyoungleaf_1", "wsyoungleaf_2", "wsyoungleaf_3", "wsroot_1", "wsroot_2", "wsroot_3", "wsberry_1", "wsberry_2", "wsberry_3")

# run DESeq on count matrices

info <-
  DataFrame(
    condition = factor(c("leaf", "leaf", "leaf", "youngleaf", "youngleaf", "youngleaf", "root", "root", "root", "berry", "berry", "berry")),
    batch = factor(c("rep1", "rep2", "rep3", "rep1", "rep2", "rep3", "rep1", "rep2", "rep3", "rep1", "rep2", "rep3")),
    row.names = colnames(counts_r)
  )

ds <-
  DESeqDataSetFromMatrix(
    countData = counts_r,
    colData = info,
    design = ~ batch + condition
  )

rld <- rlog(ds, blind = FALSE)

# choose only clusters from DESeq output and relabel

mat <- assay(rld)[c(cluster_1, cluster_2), ]
row.names(mat) <- c(labels_1, labels_2)
mat_means <- mat - rowMeans(mat)

# create heatmaps

svg(file = "clusters.svg")
ht = Heatmap(
  mat,
  name = "log2",
  cluster_rows = F,
  cluster_columns = F,
  row_split = data.frame(rep(c(
    "Cluster 1", "Cluster 2"
  ), c(19, 17))),
  column_labels = gt_render(c_labels, padding = unit(c(5, 0, 0, 0), "pt")),
  column_title = "Ws Cluster Read Counts"
)
draw(ht)
dev.off()

svg(file = "clusters_centered.svg")
ht_m = Heatmap(
  mat_means,
  name = "log2",
  cluster_rows = F,
  cluster_columns = F,
  row_split = data.frame(rep(c(
    "Cluster 1", "Cluster 2"
  ), c(19, 17))),
  column_labels = gt_render(c_labels, padding = unit(c(5, 0, 0, 0), "pt")),
  column_title = "Ws Cluster Read Counts (centered)"
)
draw(ht_m)
dev.off()

# CLEAN UP #################################################
# Clear environment
rm(list = ls())

# Clear packages
detach("package:datasets", unload = TRUE)

# Clear mind :)
