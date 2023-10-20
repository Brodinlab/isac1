library(dplyr)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(rstatix)
library(readr)
library(EnhancedVolcano)

`%nin%` = Negate(`%in%`)

# Fig7f ---------------------------------
# Zhang, L., Yu, X., Zheng, L. et al. Lineage tracking reveals dynamic relationships of T cells in colorectal cancer. Nature 564, 268–272 (2018).

dirpath <- 'Figures/fig7f'
dir.create(dirpath)

obj_Zhang <- LoadH5Seurat('Data/dataset_Zhang.h5Seurat')
Idents(obj_Zhang) <- 'Clonal.status'
res <- FindMarkers(obj_Zhang, ident.1 = 'Clonal', ident.2 = 'NoClonal',
                   logfc.threshold = 0, min.pct = 0.5)
top_p_genes <- res %>% top_n(-50, p_val_adj) %>% rownames()
interested_genes <- c('LAG3', 'HAVCR2', 'CTLA4', 'PDCD1', 'CX3CR1')
lab = unique(c(top_p_genes, interested_genes))

EnhancedVolcano(res, x='avg_log2FC', y='p_val_adj', lab = rownames(res), 
                pCutoff = 1e-02, FCcutoff = NA, drawConnectors = TRUE,
                selectLab = lab,
                title = 'memory CD8 T cells in Blood (Colorectal)',
                xlim=c(-3,3),
                col= c('grey', 'black', 'red', 'orange'),
                arrowheads=F,
                subtitle = paste0('expanded(', table(obj_Zhang[[]]$Clonal.status)['Clonal'],  
                                  ') vs non-expanded(', table(obj_Zhang[[]]$Clonal.status)['NoClonal'], ') cells')) + 
  theme_minimal()
ggsave(file.path(dirpath, 'volcano mCD8T in blood (data_Zhang) expanded vs non-expanded (Colorectal).pdf'), width = 10, height = 10)
