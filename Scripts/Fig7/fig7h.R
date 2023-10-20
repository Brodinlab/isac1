library(dplyr)
library(ggplot2)
library(rstatix)
library(readr)

# Fig7h ----------------------------------

df <- read_csv('Data/adult_isac_expansion_ratio.csv')
dirpath <- 'Figures/fig7h'
dir.create(dirpath)

ggplot(df, aes(x=factor(tumor_type, levels=unique(tumor_type)), y=fraction, color=group)) + 
  geom_bar(stat = 'summary', fun=mean, aes(fill=group)) +
  geom_jitter(color='black', width=0.2) +
  theme_bw()
ggsave(file.path(dirpath, 'Fraction of expanded cells_barplot.pdf'), width = 5, height = 4)

df %>% ungroup() %>% t_test(fraction ~ group) # added manually to the figure
