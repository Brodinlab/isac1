## figure 4

library(maftools)
library(tidyverse)
library(ggsignif)
library(data.table)
rm(list=ls())

current_time <- Sys.time()
current_time <- gsub(" ", "_", current_time)
current_time <- gsub(":", "-", current_time)

## prepare the meta table
meta_data <- read.csv("all_cancer_type_mutation_group_age_included_20230524_include_somatic_tmb_and_ISAC046_tumor_levels_2023-06-30_18-11-58.csv")

## read the maf files

folder_path <- "~/maf_vep_pass_subset/" ## files in this folder include vcf files got from sarek pipelines, then filtered by PASS
                                        ## followed by subset to Agilent Sureselect kit, then tranform to maf by vcf2maf tools
                                        ## the filtered and merged maf files can be found in data folder. Code from line 17 to line 35 just to show the filter process
                                        ## the code for the plot can run from line 38
file_names <- list.files(folder_path)
merged_data <- data.frame()
for (file_name in file_names) {
  file_path <- file.path(folder_path, file_name)
  file_data <- data.table::fread(file_path) 
  merged_data <- rbind(merged_data, file_data) 
}

classification <- unique(merged_data$Variant_Classification)

need_delete <- c("Intron","5'UTR", "3'UTR","5'Flank", "3'Flank", "IGR", "Silent", "RNA") ## use the same as maftools packages

## read filtered maf files
maf_m <- lapply(paste("~/maf_vep_pass_subset/",list.files("~/maf_vep_pass_subset/"), sep = ""), read.maf, clinicalData=meta_data, verbose = FALSE, vc_nonSyn = classification)
maf_merge <- merge_mafs(maf_m, vc_nonSyn = classification)
maf_merge <- subsetMaf(maf_merge, query = "!Variant_Classification %in% need_delete") 
maf_merge <- subsetMaf(maf_merge, query = "n_depth > 7")
write.mafSummary(maf_merge, basename = "isac1_filtered_maf")

maf_merge <- read.maf("../Data/isac1_filtered_maf_maftools.maf", verbose = F)

##---------- Fig 4a
isac.mutload = tcgaCompare(maf = new_maf, cohortName = 'ISAC', logscale = TRUE,primarySite = F, capture_size = 35.8, rm_zero = FALSE, medianCol = "#7F3F98",bg_col = c("#EDF8B1", "white"),cohortFontSize = 1, axisFontSize = 1
)


##---------- Fig 4b
tmb_exom <- tmb(maf_merge, captureSize = 35.8)
meta_table_subset <- maf_merge@clinical.data
tmb_merge <- merge(tmb_exom,meta_table_subset, by.x = "Tumor_Sample_Barcode", by.y = "Tumor_Sample_Barcode")
table_median <- tmb_merge %>% group_by(tumor_level2) %>% summarise(median = median(total_perMB)) %>% arrange(median)

tmb_merge$tumor_level2 <- factor(tmb_merge$tumor_level2, levels = table_median$tumor_level2)
tmb_merge[tmb_merge$Tumor_Sample_Barcode %in% c("810T", "P13714_104", "P26101_102", "P13714_126", "P25501_136"),]$Occassion <- "Primary_0"
ggplot(tmb_merge, mapping = aes(x=tumor_level2, y=total_perMB)) + geom_boxplot(fill = "#91928E", lwd = 0.05) +
  geom_point(position = "jitter", aes(color=Occassion),size =2) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        panel.border = element_blank(), 
        panel.grid.major  = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15, colour = "black")) +   
  scale_color_manual(values=c("#394A89", "black", "#D6A449", "#522C6E")) +
  xlab("Tumor type") + ylab("TMB (per MB)") +
  labs(color = NULL) +
  scale_y_log10()

ggsave("exom_tmb_distribution_20231018.pdf",height = 8, width = 7.2)
