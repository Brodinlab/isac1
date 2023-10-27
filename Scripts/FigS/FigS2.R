# Script to generate figure S2:

# Load packages
library(tidyverse)
library(magrittr)

# Read healthy reference data
data <- read.csv("Data/healthy_reference_major_cell_freq.csv")

# Read tumor patient protein data
protein <- read.csv("Data/isac1_baseline_protein_NPX.csv", row.names = 1, check.names = FALSE)

# Read tumor patient metadata table
meta <- read.csv2("Data/isac1_metadata.csv") 
meta %<>% filter(baseline == "yes")
meta$olink_id <- as.character(meta$olink_id)

# Check sample order consistency in protein data and metadata
identical(rownames(protein), meta$olink_id)

# Figure S2 a&b&c -----
ggplot(data, aes(x = log10(age_month), y = frequency, color = cohort)) +
  facet_wrap(~celltype, scales = "free") +
  geom_point(data = subset(cytof_plot, cohort == "Healthy"), color = "grey") +
  geom_smooth(data = subset(cytof_plot, cohort == "Healthy"), method = "loess", se = TRUE, span = 1, color = "black") +
  geom_point(data = subset(cytof_plot, cohort == "ISAC"), color = "orange") +
  labs(x = "Log10(Age in months)", y = "Fraction of cells") + 
  theme(text = element_text(color = "black"), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Figure S2 d -----
data.frame(protein = colnames(protein), cor = cor(protein, meta$age_months)) %>%
  na.omit() %>%
  ggplot(aes(x = cor, y = reorder(protein, cor))) +
  geom_segment(aes(x = 0, xend = cor, y = reorder(protein, cor), yend = reorder(protein, cor)), color = "black") +
  geom_point(aes(fill = ifelse(cor > 0, "Positive", "Negative")), color = "black", shape = 21, size = 4, show.legend = FALSE) +
  scale_fill_manual(values = c("Positive" = "black", "Negative" = "white")) +
  geom_vline(xintercept = 0, color = "black") +
  labs(x = "Correlation coefficient protein vs. age)", y = "Proteins") +
  theme(text = element_text(color = "black"), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
