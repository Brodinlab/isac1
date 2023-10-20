library(dplyr)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(rstatix)
library(readr)
library(EnhancedVolcano)

`%nin%` = Negate(`%in%`)
# Fig7g ------------------------------------
# Guo, X., Zhang, Y., Zheng, L. et al. Global characterization of T cells in non-small-cell lung cancer by single-cell sequencing. Nat Med 24, 978–985 (2018). 

dirpath <- 'Figures/fig7g'
dir.create(dirpath)

obj_Guo <- LoadH5Seurat('Data/dataset_Guo.h5Seurat')

Idents(obj_Guo) <- 'Clone.Status'
res <- FindMarkers(obj_Guo, ident.1 = 'Clonal', ident.2 = 'NoClonal',
                   logfc.threshold = 0, min.pct = 0.5)
top_p_genes <- res %>% top_n(-50, p_val_adj) %>% rownames()
interested_genes <- c('LAG3', 'HAVCR2', 'CTLA4', 'PDCD1', 'CX3CR1')
lab = unique(c(top_p_genes, interested_genes))

EnhancedVolcano(res, x='avg_log2FC', y='p_val_adj', lab = rownames(res), 
                pCutoff = 1e-02, FCcutoff = NA, drawConnectors = TRUE,
                selectLab = lab,
                title = 'memory CD8 T cells in Blood (Lung)',
                # xlim=c(-3,3),
                col= c('grey', 'black', 'red', 'orange'),
                arrowheads=F,
                subtitle = paste0('expanded(', table(obj_Guo[[]]$Clone.Status)['Clonal'],  
                                  ') vs non-expanded(', table(obj_Guo[[]]$Clone.Status)['NoClonal'], ') cells')) + 
  theme_minimal()
ggsave(file.path(dirpath, 'volcano mCD8T in blood (data_Guo) expanded vs non-expanded (Lung).pdf'), width = 10, height = 10)
