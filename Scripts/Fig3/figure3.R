# Script to generate figure 3:

# Load packages
library(tidyverse)
library(magrittr)
library(vite)
library(igraph)
library(ggraph)
library(ggrepel)

# Read metadata table
meta <- read.csv2("data/isac1_metadata.csv")
meta %<>% filter(baseline == "yes")
meta$grid_id <- as.character(meta$grid_id)
meta$olink_id <- as.character(meta$olink_id)

# Read cell cluster frequency data
cell <- read.csv("data/isac1_baseline_cell_cluster_freq.csv", row.names = 1, check.names = FALSE)

# Read cell cluster median marker expression data
marker <- read.table("data/isac1_baseline_cell_cluster_median_marker_expression.txt", sep = "\t", header = TRUE)

# Read protein expression data
protein <- read.csv("data/isac1_baseline_protein_NPX.csv", row.names = 1)

# Check sample order consistency in data and metadata
identical(rownames(cell), meta$grid_id)
identical(rownames(protein), meta$olink_id)

# Figure 3a & 3c -----

# For figure 3a: filter neuroblastoma patients
data <- cell %>%
  cbind(select(meta, tumor_level1, metastasis)) %>%
  filter(tumor_level1 == "Neuroblastoma") %>%
  mutate(group = metastasis) %>%
  select(-tumor_level1, -metastasis) %>%
  select_if(~!all(is.na(.)))

# For figure 3c: filter patients with medication
data <- cell %>%
  cbind(select(meta, medication, neutropenic_fever)) %>%
  filter(medication == "Yes") %>%
  mutate(group = neutropenic_fever) %>%
  select(-medication, -neutropenic_fever) %>%
  select_if(~!all(is.na(.)))

# Separate two comparison groups
group1 <- data %>% filter(group == "Yes") %>% select(-group)
group2 <- data %>% filter(group == "No") %>% select(-group)

# Calculate the log2 fold change 
log2fc <- log2(colMeans(group1, na.rm = TRUE) / colMeans(group2, na.rm = TRUE))

# Perform Wilcoxon test for each feature
pvalue <- sapply(1:ncol(group1), function(i) wilcox.test(group1[, i], group2[, i])$p.value)

# Combine results
df <- data.frame(
  feature = colnames(group1),
  log2fc = log2fc,
  pvalue = pvalue
)

# Build force-directed graph from cell cluster median marker expression
marker_df <- marker %>% mutate(feature = celltype) %>% left_join(df)
set.seed(824)
G <- build_graph(marker_df, col.names = colnames(marker_df)[1:32], filtering_T = 5)
for(i in names(marker_df)){
  G <- set.vertex.attribute(G, name = i, value = marker_df[, i])
}
cc <- multilevel.community(G)
V(G)$community_id <- as.character(cc$membership) 
V(G)$type <- "cluster"
V(G)$Label <- paste("c", V(G)$cellType, sep = "")

# Visualize force-directed graph
set.seed(824)
l_g <- create_layout(G, layout="fr")
l_g$x <- V(G)$x
l_g$y <- V(G)$y
ggraph(l_g) +
  geom_edge_link(alpha = .1, check_overlap = F) + 
  geom_node_point(aes(size = percentage), fill = "grey", shape = 21,  alpha = 1) +
  geom_node_point(data = subset(l_g, pvalue < 0.1), 
                  aes(size = percentage, fill = ifelse(log2fc > 0, "Positive", "Negative")), 
                  shape = 21,  alpha = 1) +
  scale_size_continuous(range = c(5, 20)) +
  scale_fill_manual(values = c("Positive" = "#65A479FF", "Negative" = "#D3BA68FF"), name = "Log2 Fold Change") +
  geom_node_text(data = subset(l_g, pvalue < 0.1),
                 aes(label = cluster), size = 4, color = 'white', show.legend = FALSE, repel = F) +
  theme_graph(base_family = 'Helvetica') 

#  Figure 3b & 3c -----

# For figure 3b: filter neuroblastoma patients
data <- protein %>%
  cbind(select(meta, tumor_level1, metastasis)) %>%
  filter(tumor_level1 == "Neuroblastoma") %>%
  mutate(group = metastasis) %>%
  select(-tumor_level1, -metastasis) %>%
  select_if(~!all(is.na(.)))

# For figure 3c: filter patients with medication
data <- protein %>%
  cbind(select(meta, medication, neutropenic_fever)) %>%
  filter(medication == "Yes") %>%
  mutate(group = neutropenic_fever) %>%
  select(-medication, -neutropenic_fever) %>%
  select_if(~!all(is.na(.)))

# Separate two comparison groups
group1 <- data %>% filter(group == "Yes") %>% select(-group)
group2 <- data %>% filter(group == "No") %>% select(-group)

# Calculate the log2 fold change 
log2fc <- log2(colMeans(group1, na.rm = TRUE) / colMeans(group2, na.rm = TRUE))

# Perform Wilcoxon test for each feature
pvalue <- sapply(1:ncol(group1), function(i) wilcox.test(group1[, i], group2[, i])$p.value)

# Combine results
df <- data.frame(
  feature = colnames(group1),
  log2fc = log2fc,
  pvalue = pvalue
)

# Volcano plot
ggplot(df, aes(x = log2fc, y = -log10(pvalue))) +
  geom_point(aes(color = ifelse(log2fc > 0, "Positive", "Negative")), size = 2) +
  geom_text_repel(data = subset(df, pvalue < 0.5), aes(label = feature)) + 
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  scale_color_manual(values = c("Positive" = "#65A479FF", "Negative" = "#D3BA68FF"), name = "Log2 Fold Change") +
  labs(x = "Log2 Fold Change", y = "-log10(P-Value)") +
  theme(panel.background = element_blank())
