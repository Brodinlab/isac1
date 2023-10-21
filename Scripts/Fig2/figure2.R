# Script to generate figure 2:

# Load packages
library(tidyverse)
library(magrittr)
library(robCompositions)
library(vegan)
library(pheatmap)

# Read metadata table
meta <- read.csv2("data/isac1_metadata.csv")
meta %<>% filter(baseline == "yes")
meta$grid_id <- as.character(meta$grid_id)
meta$olink_id <- as.character(meta$olink_id)

# Read cell composition data
cell <- read.csv("data/isac1_baseline_cell_subtype_freq.csv", row.names = 1, check.names = FALSE)

# Read protein expression data
protein <- read.csv("data/isac1_baseline_protein_NPX.csv", row.names = 1, check.names = FALSE)

# Check sample order consistency in data and metadata
identical(rownames(cell), meta$grid_id)
identical(rownames(protein), meta$olink_id)

# Figure 2b -----

# Compute the Aitchison distances between samples
aitchison_dist <- aDist(cell)

# Classical multidimensional scaling
mds <- cmdscale(aitchison_dist)

# Add metadata
mds_df <- cbind(mds, meta)

# MDS plot
ggplot(mds_df, aes(x = `1`, y = `2`, color = age_months)) +
  geom_point(size = 3) +
  paletteer::scale_color_paletteer_c("ggthemes::Purple") +
  labs(x = "MDS1", y = "MDS2", title = "MDS of immune cell composition") +
  theme(panel.background = element_blank()) 

# Figure 2c -----

# Remove proteins with variance lower than 5% of the total variance
protein_variance <- apply(protein, 2, var)
which(protein_variance < 0.05*mean(var(protein, na.rm = TRUE), na.rm = TRUE))

# Principle component analysis
pca <- prcomp(protein[,colSums(is.na(protein)) == 0], scale. = TRUE) # Remove proteins with any NAs

# PCA plot
autoplot(pca, data = meta, color = "age_months", size = 3) +
  paletteer::scale_color_paletteer_c("ggthemes::Purple") +
  theme(panel.background = element_blank())

# Figure 2d & 2e -----

# For figure 2d:
data <- cell
batch <- meta$cytof_batch
vegdist_method <- "robust.aitchison"

# For figure 2e:
data <- protein
batch <- meta$olink_batch
vegdist_method <- "euclidean"

# Replace NAs with zeros
data[is.na(data)] <- 0

# Perform PERMANOVA analysis on each experimental batch
permanova <- list()

# Loop over unique batches
for (batch_i in  unique(batch)) {
  batch_row <- which(batch == batch_i)
  # Skip batches encountering errors
  try({
    permanova_result <- adonis2(data[batch_row,] ~ sex + age_months + tumor_level1 + metastasis,
                                meta[batch_row,], method = vegdist_method, by = "margin")
    if (nrow(permanova_result) == 3) {
      permanova[[batch_i]] <- NA
    } else {
      permanova[[batch_i]] <- data.frame(batch = batch_i,
                                         factor = c("Sex", "Age", "Tumor type", "Metastasis", "Residual", "Total"),
                                         R2 = permanova_result$R2,
                                         pvalue = permanova_result$`Pr(>F)`)
    }
    
  }, silent = TRUE)
  
  if (is.na(permanova[batch_i])) {
    permanova[[batch_i]] <- NA
  }
}

# Bind results together
permanova <- do.call(rbind, permanova) %>%
  group_by(batch) %>%
  mutate(n_NA = sum(is.na(pvalue))) %>%
  # Remove batches with complete enumeration
  filter(n_NA != length(unique(batch))) %>%
  # Remove total 
  filter(factor != "Total") %>%
  # Recalculate R2 to sum up to 1
  mutate(sum_R2 = sum(R2)) %>%
  mutate(R2 = R2/sum_R2) %>%
  select(-n_NA, -sum_R2) 

# Set factor levels
permanova$factor <- factor(permanova$factor, levels = c("Sex", "Metastasis", "Age", "Tumor type", "Residual"))

# Visualization
permanova %>%
  ggplot(aes(batch, R2, fill = factor)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  labs(x = "Batch", y = "Percentage of variance", 
       title = "Variance explained by each factor in each batch") +
  scale_fill_manual(values = paletteer::paletteer_d("fishualize::Chaetodon_sedentarius", direction = -1)) +
  theme(panel.background = element_blank()) 

# Figure 2f & 2g -----

# For figure 2f:
data <- cell
dist_fun <- function(data_median){aDist(data_median)}

# For figure 2g:
data <- protein
dist_fun <- function(data_median){dist(data_median, method = "euclidean")}

# Exclude tumor types with 2 or less patients
tumor_counts <- table(meta$tumor_level1)
tumor_exclude <- names(tumor_counts[tumor_counts %in% c(1,2)])

# Compute median value of features per tumor type 
data_median <- data %>% 
  mutate(tumor_level1 = meta$tumor_level1) %>%
  filter(!tumor_level1 %in% tumor_exclude) %>%
  gather(feature, value, -tumor_level1) %>%
  group_by(tumor_level1, feature) %>%
  mutate(median = median(value)) %>%
  select(-value) %>%
  distinct() %>%
  spread(feature, median) %>%
  column_to_rownames("tumor_level1")

# Compute distances between tumor types
dist <- dist_fun(data_median) %>% as.matrix()

# Reshape dist table
dist_df <- dist %>%
  as.table() %>%
  as.data.frame()
colnames(dist_df) <- c("From", "To", "Distance")

# Plot heatmap with hierarchical tree
heatmap <- pheatmap(dist, cluster_rows = TRUE, cluster_cols = TRUE,
                    angle_col = 315,
                    fontsize = 15,
                    display_numbers = TRUE)

# Reorder tumor types according to hierarchical tree
dist_df$From = factor(dist_df$From, levels = heatmap$tree_row$labels[heatmap$tree_row$order])
dist_df$To = factor(dist_df$To, levels = rev(heatmap$tree_row$labels[heatmap$tree_row$order]))

# Plot dotplot
dotplot <- ggplot(dist_df, aes(x = From, y = To, color = -Distance, size = -Distance)) +
  geom_point() +
  scale_color_gradient(low = "white", high = "gray30") +
  scale_size_continuous(range = c(1, 15)) +
  scale_y_discrete(position = "right") +
  labs(x = "", y = "", title = "Pairwise distances between tumor types") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        plot.title = element_text(size = 10),
        panel.background = element_blank())
combined_legend <- guides(
  color = guide_legend(title = "Distance", override.aes = list(shape = 16)),
  size = guide_legend(title = "Distance", override.aes = list(shape = 16))
)
dotplot + combined_legend
