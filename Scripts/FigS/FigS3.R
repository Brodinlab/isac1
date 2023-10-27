library(dplyr)
library(readr)
source('Scripts/func/GLIPH_related_plots.R')
source("Scripts/func/clone_expansion_plots.R")


data_scTCR <- read_csv('Data/scTCR_baseline.csv.gz')
data_TRA <- read_csv('Data/GLIPH2_TRA_annotated.csv.gz')
data_TRB <- read_csv('Data/GLIPH2_TRB_annotated.csv.gz')
clone_id_map <- read_csv('Data/clone_id_map.csv.gz')
clone_exp <- read_csv(file = 'Data/clone_expansion.csv.gz') %>%
  mutate(clone_id = as.character(clone_id),
         patient_id = sub('_.*$', '', Sample_Name))

nt2aa <- data_scTCR %>% select(CDR3_concat, CDR3aa_concat) %>% distinct()

plot.gliph(data_TRA %>% filter(condition %in% c('1', 'v1')), 
           data_TRB %>% filter(condition %in% c('1', 'v1')), 
           clone_exp %>% filter(grepl('(_v1$)|(_1$)', Sample_Name)), 
           nt2aa, 
           clone_thre=1) 
# this can be a bit different from the one in the figure because I didn't fix the random seed.
