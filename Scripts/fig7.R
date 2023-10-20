library(dplyr)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(rstatix)
library(readr)
library(EnhancedVolcano)

`%nin%` = Negate(`%in%`)

# Fig7a-d ------------------------------------
obj <- LoadH5Seurat('Data/seurat_results_CD8T.h5Seurat')
expand_sc <- read_csv('Data/expand_sc.csv.gz')


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

for (cur_tumor_type in unique(obj[[]]$tumor_type)){
  # cur_tumor_type = 'Wilms'
  obj_sub <- subset(obj, subset = tumor_type == cur_tumor_type)
  volcano_from_obj(obj_sub, expand_sc, title = paste0('CD8T of ', cur_tumor_type))
  ggsave(paste0('Figures/fig7/CD8T_', cur_tumor_type ,'_expansion.pdf'), width = 15, height = 15)
}

# Fig7e ------------------------------------
# Wu, T.D., Madireddi, S., de Almeida, P.E. et al. Peripheral T cell expansion predicts tumour infiltration and clinical response. Nature 579, 274–278 (2020).

df_Wu <- read_csv('Data/blood mCD8T from dataset_Wu.csv.gz')
panel_Wu <- colnames(df_Wu)[1:3734]

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
ggsave('Figures/fig7/volcano mCD8T in blood (data_Wu) expanded vs non-expanded (all).pdf', width = 10, height = 10)

# Fig7f ---------------------------------
# Zhang, L., Yu, X., Zheng, L. et al. Lineage tracking reveals dynamic relationships of T cells in colorectal cancer. Nature 564, 268–272 (2018).

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
ggsave('Figures/fig7/volcano mCD8T in blood (data_Zhang) expanded vs non-expanded (Colorectal).pdf', width = 10, height = 10)

# Fig7g ------------------------------------
# Guo, X., Zhang, Y., Zheng, L. et al. Global characterization of T cells in non-small-cell lung cancer by single-cell sequencing. Nat Med 24, 978–985 (2018). 

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
ggsave('Figures/fig7/volcano mCD8T in blood (data_Guo) expanded vs non-expanded (Lung).pdf', width = 10, height = 10)

# Fig7h ----------------------------------

df <- read_csv('Data/adult_isac_expansion_ratio.csv')

ggplot(df, aes(x=factor(tumor_type, levels=unique(tumor_type)), y=fraction, color=group)) + 
  geom_bar(stat = 'summary', fun=mean, aes(fill=group)) +
  geom_jitter(color='black', width=0.2) +
  theme_bw()
ggsave('Figures/fig7/Fraction of expanded cells_barplot.pdf', width = 5, height = 4)

df %>% ungroup() %>% t_test(fraction ~ group) # added manually to the figure


