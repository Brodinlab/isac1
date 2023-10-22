# Script to generate figure S1:

# Load packages
library(pheatmap)

# Read cell cluster median marker expression data
marker <- read.table("Data/isac1_baseline_cell_cluster_median_marker_expression.txt", sep = "\t", header = TRUE)

# Plot heatmap
pheatmap(marker[,1:32],
         color = paletteer::paletteer_d("rcartocolor::PurpOr", n = 100, type = 'continuous'), 
         scale = "column",
         labels_row = marker$celltype,
         display_numbers = FALSE,
         angle_col = 45,
         breaks = seq(0, 2, 2/90),
         main = "Median marker expression per cluster (z-score transformed)")
