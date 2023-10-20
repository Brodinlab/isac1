library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(EnhancedVolcano)
library(Seurat)
library(SeuratDisk)

source("Scripts/func/clone_expansion_plots.R")

clone_id_map <- read_csv('Data/clone_id_map.csv.gz')
virus_specific_list <- read_csv(file = 'Data/virus_specific_clone_id.csv') # obtain from GLIPH2 results
data_scTCR <- read_csv('Data/scTCR_data_ISAC_02.csv.gz') %>%
  mutate(seurat_clusters = as.character(seurat_clusters))
obj <- LoadH5Seurat('Data/seurat_results_fig6.h5Seurat')

fig_dir <- 'Figures/fig6'
patient <- 'ISAC02'


# fig 6a-b --------------
for (sub_name in c('CD4T', 'CD8T', 'gdT')){
  clone_exp_sub <- get_clone_exp_sub(data_scTCR, sub_name) %>%
    inner_join(clone_id_map, by='CDR3_concat')
  g_list1 <- lapply(unique(data_scTCR$Sample_Name),function(x){clone_expansion_donut(x, clone_exp_sub)})
  ggarrange(plotlist = g_list1, ncol = 3, nrow = 5) %>%
    ggexport(filename = file.path(fig_dir, paste0('clone_expansion_', sub_name, '.pdf')), 
             width = 10, height = 20)
  g <- clone_expansion_alluvium(patient, clone_exp_sub) + labs(title = paste0(patient, '_', sub_name))
  ggsave(plot = g, filename = file.path(fig_dir, paste0('top_clone_changes_', sub_name, '.pdf')))
}

# fig 6c --------------

sub_name <- 'CD8T'
clone_exp_sub <- get_clone_exp_sub(data_scTCR, sub_name) %>%
  inner_join(clone_id_map, by='CDR3_concat')
g <- clone_expansion_alluvium(patient, clone_exp_sub, clone_label = T) + labs(title = paste0(patient, '_', sub_name)) + guides(fill='none', color='none')

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

g3 <- ggplot(df, aes(x = time_point, stratum = clone_id, alluvium = clone_id, 
                     #y= clone_ratio, 
                     y = relative_clone_ratio,
                     fill = trend, color=clone_id)) +
  geom_alluvium() +
  geom_stratum(size = 0.1) +
  guides(color = "none") + # hide the legend of 'color'
  labs(title = patient)

# fig 6d --------------

tmp <- obj[[]] %>% left_join(df %>% select(clone_id, trend) %>% distinct(), by='clone_id') %>% select(trend)
rownames(tmp) <- obj[[]]$unique_index
obj$trend <- tmp
Idents(obj) <- 'trend'
res <- FindMarkers(obj, ident.1 = 'increase', ident.2 = 'decrease',
                   logfc.threshold = 0)

g4 <- EnhancedVolcano(res, x='avg_log2FC', y='p_val', lab = rownames(res), 
                      pCutoff = 1e-04, FCcutoff = 1, drawConnectors = TRUE,
                      selectLab = res %>% filter(abs(avg_log2FC) > 2 | p_val < 1e-04) %>% rownames(),
                      title = paste0(sub_name, ' of ', patient),
                      subtitle = 'increase vs decrease')

ggarrange(g, g3, g4, ncol = 3, common.legend=T) %>% ggexport(filename = paste0(fig_dir, '/top_clone_changes_CD8T_trend.pdf'), width = 20, height = 8)
