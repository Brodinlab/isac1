library(dplyr)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(rstatix)
library(readr)
library(EnhancedVolcano)

obj <- LoadH5Seurat('Data/seurat_results_CD8T.h5Seurat')
expand_sc <- read_csv('Data/expand_sc.csv.gz')
dirpath <- 'Figures/fig7b'
dir.create(dirpath)

volcano_from_obj <- function(oj, exp_sc, title){
  
  tmp <- oj[[]] %>% left_join(exp_sc, by = 'unique_index') %>%
    tidyr::replace_na(list(expansion='non-expanded'))
  tmp2 <- tmp$expansion
  names(tmp2) <- tmp$unique_index
  oj$expansion <- tmp2
  Idents(oj) <- 'expansion'
  res <- FindMarkers(oj, ident.1 = 'expanded', ident.2 = 'non-expanded',
                     logfc.threshold = 0, min.pct = 0)
  top_p_genes <- res %>% top_n(-50, p_val_adj) %>% rownames()
  interested_genes <- c('LAG3', 'HAVCR2', 'CTLA4', 'PDCD1')
  lab = unique(c(top_p_genes, interested_genes))
  
  EnhancedVolcano(res, x='avg_log2FC', y='p_val_adj', lab = rownames(res), 
                  pCutoff = 1e-02, FCcutoff = NA, drawConnectors = TRUE,
                  selectLab = lab,
                  title = title,
                  xlim=c(-1.5,1.5),
                  col= c('grey', 'black', 'orange', 'red'),
                  arrowheads=F,
                  subtitle = paste0('expanded vs non-expanded cells')) + 
    theme_minimal()
}


cur_tumor_type = 'High mutation Neuroblastoma'
obj_sub <- subset(obj, subset = tumor_type == cur_tumor_type)
volcano_from_obj(obj_sub, expand_sc, title = paste0('CD8T of ', cur_tumor_type))
ggsave(paste0(dirpath,'/CD8T_', cur_tumor_type ,'_expansion.pdf'), width = 15, height = 15)

