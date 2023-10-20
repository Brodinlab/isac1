library(dplyr)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(rstatix)
library(readr)
library(EnhancedVolcano)

`%nin%` = Negate(`%in%`)

# Fig7e ------------------------------------
# Wu, T.D., Madireddi, S., de Almeida, P.E. et al. Peripheral T cell expansion predicts tumour infiltration and clinical response. Nature 579, 274–278 (2020).

df_Wu <- read_csv('Data/blood mCD8T from dataset_Wu.csv.gz')
panel_Wu <- colnames(df_Wu)[1:3734]
dirpath <- 'Figures/fig7e'
dir.create(dirpath)

volcano_tmp <- function(df){
  tmp <- df %>% group_by(expansion) %>% summarise(across(all_of(panel_Wu), mean))
  od <- tmp$expansion
  tmp <- t(tmp[, panel_Wu])
  colnames(tmp) <- od
  resFC <- tmp %>% as.data.frame() %>% mutate(logFC = expanded - non_expanded)
  resFC$gene <- rownames(resFC)
  tmp <- df %>% tidyr::pivot_longer(values_to = 'level', names_to = 'gene', cols = all_of(panel_Wu))
  resP <- tmp %>% 
    group_by(gene) %>% 
    t_test(level ~ expansion) %>% 
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj")
  res <- left_join(resFC, resP, by = 'gene')
  top_p_genes <- (res %>% top_n(-50, p.adj))$gene
  interested_genes <- c('LAG3', 'HAVCR2', 'CTLA4', 'PDCD1', 'CX3CR1')
  lab = unique(c(top_p_genes, interested_genes))
  
  EnhancedVolcano(res, x='logFC', y='p.adj', lab = res$gene, 
                  pCutoff = 1e-02, FCcutoff = NA, drawConnectors = TRUE, arrowheads = F, 
                  selectLab = lab,
                  xlim=c(-1,1),
                  col= c('grey', 'black', 'red', 'orange'),
                  subtitle = paste0('expanded(', table(df$expansion)['expanded'],  
                                    ') vs non-expanded(', table(df$expansion)['non_expanded'], ') cells')) + 
    theme_minimal()
}

volcano_tmp(df_Wu) + ggtitle('memory CD8 T cells in Blood (all)')
ggsave(file.path(dirpath, 'volcano mCD8T in blood (data_Wu) expanded vs non-expanded (all).pdf'), width = 10, height = 10)
