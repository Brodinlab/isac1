# Script to generate figure 6e:

# Load packages
library(tidyverse)
library(ggrepel)

# Read protein data
protein <- read.csv("Data/isac002_protein.csv", row.names = 1, check.names = F)

# Calculate log2 fold change
df <- protein %>%
  rownames_to_column("protein") %>%
  mutate(log2fc = log2(`9.month`/`5.month`)) 

# Figure 6e -----
ggplot(df, aes(x = log2fc, y = reorder(protein, log2fc))) +
  geom_vline(linetype = "longdash", xintercept = 0, color = "grey") +
  geom_point(color = "black", size = 4) +
  geom_text_repel(data = subset(df, abs(log2fc) > 0.6), aes(label = protein), box.padding = 0.5) +
  labs(x = "Log2fc(9 month/5 month)", y = "Proteins") +
  theme(text = element_text(color = "black"), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
