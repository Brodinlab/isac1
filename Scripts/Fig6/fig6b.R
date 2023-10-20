library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)

source("Scripts/func/clone_expansion_plots.R")

clone_id_map <- read_csv('Data/clone_id_map.csv.gz')
data_scTCR <- read_csv('Data/scTCR_data_ISAC_02.csv.gz') %>%
  mutate(seurat_clusters = as.character(seurat_clusters))

fig_dir <- 'Figures/fig6b'
patient <- 'ISAC02'
dir.create(fig_dir)
sub_name <- 'CD8T'

clone_exp_sub <- get_clone_exp_sub(data_scTCR, sub_name) %>%
  inner_join(clone_id_map, by='CDR3_concat')
g <- clone_expansion_alluvium(patient, clone_exp_sub) + labs(title = paste0(patient, '_', sub_name))
ggsave(plot = g, filename = file.path(fig_dir, paste0('top_clone_changes_', sub_name, '.pdf')))

