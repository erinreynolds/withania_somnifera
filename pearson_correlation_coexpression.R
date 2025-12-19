#written by Jonathan Huang

library(corrplot)
library(gridGraphics)  # To convert base R plots to grid graphics
library(grid)

#import normalized counts and transpose 
normalized_counts_final <- read.table("normalized_counts_final.txt", row.names = 1, header = T)
norm_counts <- t(normalized_counts_final)

#calculate correlation matrix between every variable in norm_counts; took several minutes to run
PCC_matrix <- cor(norm_counts, method = 'pearson')

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

#make a matrix of the clusters only
clus1 <- PCC_matrix[cluster_1, cluster_1]
row.names(clus1) <- c(labels_1)
colnames(clus1) <- c(labels_1)
clus2 <- PCC_matrix[cluster_2, cluster_2]
row.names(clus2) <- c(labels_2)
colnames(clus2) <- c(labels_2)

# Set the font to italics and Arial --> the italics is not working
svg(file="cluster1_cor_matrix_updated_JH.svg")

#create correlation plot and save as svg
# tl.font = 3 only works for generating the image, but this is for some reason
# incompatible with exporting to svg
# so turn off the default thing, do some weird grid thing
corrplot(clus1, tl.col = "black", tl.pos = "n", mar = c(1, 5, 2, 1))

# Convert base R plot to grid object
grid.echo()

# Italicize column and row labels using grid.text()
## feel free to fine tune the coordinate values
## if you don't want the rotation set rot = 0
n <- ncol(clus1)
grid.text(rownames(clus1), x = unit(0.15, "npc"), y = seq(0.815, 0.16, length = n),
          just = "right", gp = gpar(fontface = "italic", cex = 1), rot = 0)

grid.text(colnames(clus1), x = seq(0.175, 0.835, length = n), y = unit(0.85, "npc"),
          just = c("left"), gp = gpar(fontface = "italic", cex = 1), rot = 45)

dev.off()

svg(file="cluster2_cor_matrix_updated_JH.svg")
# same as above=
corrplot(clus2, tl.col = "black", tl.pos = "n", mar = c(1, 5, 2, 1))
grid.echo()
n <- ncol(clus2)
grid.text(rownames(clus2), x = unit(0.15, "npc"), y = seq(0.815, 0.16, length = n),
          just = "right", gp = gpar(fontface = "italic", cex = 1), rot = 0)
grid.text(colnames(clus2), x = seq(0.175, 0.835, length = n), y = unit(0.85, "npc"),
          just = c("left"), gp = gpar(fontface = "italic", cex = 1), rot = 45)
dev.off()

# CLEAN UP #################################################
# Clear environment
#rm(list = ls())

# Clear mind :)