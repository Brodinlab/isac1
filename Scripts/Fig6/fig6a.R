library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)

source("Scripts/func/clone_expansion_plots.R")

clone_id_map <- read_csv('Data/clone_id_map.csv.gz')
data_scTCR <- read_csv('Data/scTCR_data_ISAC_02.csv.gz') %>%
  mutate(seurat_clusters = as.character(seurat_clusters))

fig_dir <- 'Figures/fig6a'
patient <- 'ISAC02'

dir.create(fig_dir)

sub_name = 'CD8T'
clone_exp_sub <- get_clone_exp_sub(data_scTCR, sub_name) %>%
  inner_join(clone_id_map, by='CDR3_concat')
g_list1 <- lapply(unique(data_scTCR$Sample_Name),function(x){clone_expansion_donut(x, clone_exp_sub)})
ggarrange(plotlist = g_list1, ncol = 3, nrow = 5) %>%
  ggexport(filename = file.path(fig_dir, paste0('clone_expansion_', sub_name, '.pdf')), 
           width = 10, height = 20)

