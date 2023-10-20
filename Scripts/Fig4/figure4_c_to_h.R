# Script to generate figure 4c to 4h:

# Load packages
library(tidyverse)
library(magrittr)
library(ggpubr)

# Read metadata table
meta <- read.csv2("data/isac1_metadata.csv")
meta$grid_id <- as.character(meta$grid_id)
meta$olink_id <- as.character(meta$olink_id)

# Filter neuroblastoma patients
meta %<>% filter(tumor_level1 == "Neuroblastoma" & baseline == "yes" & is.na(PID) == FALSE)

# Read cell composition data
cell <- read.csv("data/isac1_baseline_cell_subtype_freq.csv", row.names = 1, check.names = FALSE)
cell <- cell[meta$grid_id,]  %>% select_if(~!all(is.na(.)))

# Read protein expression data
protein <- read.csv("data/isac1_baseline_protein_NPX.csv", row.names = 1)
protein <- protein[meta$olink_id,]  %>% select_if(~!all(is.na(.)))

# Check sample order consistency in data and metadata
identical(rownames(cell), meta$grid_id)
identical(rownames(protein), meta$olink_id)

# Figure 4c & 4d & 4e -----
cell %>% 
  cbind(meta %>% select(study_id, tmb_group)) %>%
  gather(feature, value, -study_id, -tmb_group) %>%
  ggplot(aes(x = tmb_group, y = value, fill = tmb_group)) +
  facet_wrap(~feature, scales = "free") +
  stat_compare_means(method = "wilcox.test", vjust = 1) +
  geom_boxplot(alpha = 0.6) +
  geom_point() +
  geom_text(aes(label = study_id)) +
  scale_fill_manual(values = c("#D5695DFF", "#5D8CA8FF")) +
  theme(panel.background = element_blank())
  
# Figure 4f -----
protein %>%
  cbind(meta %>% select(tmb_group)) %>%
  gather(feature, value, -tmb_group) %>%
  na.omit() %>%
  group_by(feature, tmb_group) %>%
  mutate(median = median(value)) %>%
  select(-value) %>%
  distinct() %>%
  spread(tmb_group, median) %>%
  na.omit() %>%
  mutate(log2fc = log2(high/low)) %>%
  ggplot(aes(x = log2fc, y = reorder(feature, log2fc))) +
  geom_segment(aes(x = 0, xend = log2fc, y = reorder(feature, log2fc), yend = reorder(feature, log2fc)), color = "grey") +
  geom_vline(xintercept = 0, color = "grey") +
  geom_point(aes(color = ifelse(log2fc > 0, "Positive", "Negative")), size = 4, shape = 21, fill = "white", show.legend = FALSE) +
  scale_color_manual(values = c("Positive" = "brown", "Negative" = "#D3BA68FF")) +
  labs(x = "Log2fc(high mutation/low mutation)", y = "Proteins") +
  theme(panel.background = element_blank())

# Figure 4g -----
cell %>%
  cbind(meta %>% select(ros)) %>%
  gather(feature, value, -ros) %>%
  na.omit() %>%
  group_by(feature, ros) %>%
  mutate(median = median(value)) %>%
  select(-value) %>%
  distinct() %>%
  spread(ros, median) %>%
  na.omit() %>%
  mutate(log2fc = log2(Yes/No)) %>%
  ggplot(aes(x = log2fc, y = reorder(feature, log2fc))) +
  geom_segment(aes(x = 0, xend = log2fc, y = reorder(feature, log2fc), yend = reorder(feature, log2fc)), color = "grey") +
  geom_vline(xintercept = 0, color = "grey") +
  geom_point(aes(color = ifelse(log2fc > 0, "Positive", "Negative")), size = 4, shape = 21, fill = "white", show.legend = FALSE) +
  scale_color_manual(values = c("Positive" = "brown", "Negative" = "#D3BA68FF")) +
  labs(x = "Log2fc(ROS/non-ROS)", y = "Proteins") +
  theme(panel.background = element_blank())

# Figure 4h -----

# Remove proteins with any NA 
protein <- protein[, complete.cases(t(protein))]

# Remove proteins with low variance
mean(var(protein))
protein_variance <- apply(protein, 2, var)
which(protein_variance < 0.05*mean(var(protein)))
protein <- protein[,!protein_variance < 0.05*mean(var(protein))]

# Perform PCA
pca_res <- prcomp(protein, scale. = TRUE)

# Plot PCA result
autoplot(pca_res, data = meta, color = "ros", size = 5, loadings = TRUE, loadings.color = "lightgrey",
         loadings.label = 1, loadings.label.size = 2, loadings.label.color = "#0C2340FF", loadings.label.repel = TRUE) +
  scale_color_manual(values = c("#FFA400FF","#2C5234FF")) +
  theme(panel.background = element_blank()) +
  labs(title = "PCA of Olink data of Neuroblastoma patients: ROS vs. non-ROS")
