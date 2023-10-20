library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)

source("Scripts/func/clone_expansion_plots.R")

clone_id_map <- read_csv('Data/clone_id_map.csv.gz')
virus_specific_list <- read_csv(file = 'Data/virus_specific_clone_id.csv') # obtain from GLIPH2 results
data_scTCR <- read_csv('Data/scTCR_data_ISAC_02.csv.gz') %>%
  mutate(seurat_clusters = as.character(seurat_clusters))

fig_dir <- 'Figures/fig6c'
patient <- 'ISAC02'

dir.create(fig_dir)

sub_name <- 'CD8T'
clone_exp_sub <- get_clone_exp_sub(data_scTCR, sub_name) %>%
  inner_join(clone_id_map, by='CDR3_concat')

g2 <- trend_determination_plot(patient, clone_exp_sub, top_mod='union')

# determine manually according to g2
increase_clone <- c(48795, 37805, 2890, 45479, 53017, 33890)
decrease_clone <- c(7918, 30260, 51808, 32056, 4293, 41385, 18961, 50670, 41552, 
                    47305, 20424, 18818, 47227, 32052, 43590, 32214, 55081)


# flow by trend
df <- clone_exp_alluvium_preparation(patient, clone_exp_sub, top_mod='union')
df <- df %>% mutate(
  trend = case_when(
    clone_id %in% virus_specific_list$clone_id ~ 'virus_specific',
    clone_id %in% increase_clone ~ 'increase',
    clone_id %in% decrease_clone ~ 'decrease',
    TRUE ~ 'other'
  ))

ggplot(df, aes(x = time_point, stratum = clone_id, alluvium = clone_id, 
                     #y= clone_ratio, 
                     y = relative_clone_ratio,
                     fill = trend, color=clone_id)) +
  geom_alluvium() +
  geom_stratum(size = 0.1) +
  guides(color = "none") + # hide the legend of 'color'
  labs(title = patient)
ggsave(file.path(fig_dir, 'fig6c.pdf'))
