library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)

source("Scripts/func/clone_expansion_plots.R")

# load data
clone_id_map <- read_csv('data/clone_id_map.csv.gz')
data_scTCR <- read_csv('data/scTCR_baseline.csv.gz') %>%
  mutate(seurat_clusters = as.character(seurat_clusters))

# donut chart 
fig_dir = 'Figures/fig5c'
sample_names <- c('ISAC77_v1')

for (sub_name in c('CD4T', 'CD8T', 'gdT')){
  clone_exp_sub <- get_clone_exp_sub(data_scTCR, sub_name) %>%
    inner_join(clone_id_map, by='CDR3_concat')
  g_list1 <- lapply(sample_names,function(x){clone_expansion_donut(x, clone_exp_sub)})
  ggarrange(plotlist = g_list1, ncol = 3, nrow = 5) %>%
    ggexport(filename = file.path(fig_dir, paste0('clone_expansion_', sub_name, '.pdf')), 
             width = 10, height = 20)
}

