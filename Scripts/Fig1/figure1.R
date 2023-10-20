# Script to generate figure1:

# Load packages
library(tidyverse)
library(webr)

# Read functions
source("scripts/PieDonutCustom.R")

# Read metadata table
meta <- read.csv2("data/isac1_metadata.csv") %>% filter(is.na(subject_id) == FALSE)

# Reorder tumor type levels
meta$tumor_type_grouped <- factor(meta$tumor_type_grouped, 
                                  levels = c("Lymphoma", "Brain tumor", "Neuroblastoma", 
                                             "Retinoblastoma", "Kidney tumor", "Liver tumor", 
                                             "Bone tumor", "Soft tissue sarcoma", 
                                             "Germ cell tumors", "Others"))

# Figure 1a -----
meta %>%
  count(sex, age_year) %>%
  ggplot(aes(age_year, n, fill = sex)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  coord_flip() +
  scale_x_reverse(breaks = seq(0, 18, by = 1)) +
  scale_y_continuous(breaks = seq(0, 15, by = 1), position = "right") +
  labs(x = "Number of patients", y = "Age (Years)") +
  scale_fill_manual(values = c("darkgrey","#A7D3D4FF")) +
  theme(panel.background = element_blank())

# Figure 1b -----
meta %>% 
  group_by(tumor_type_grouped, tumor_level2) %>%
  tally() %>%
  PieDonutCustom(aes(tumor_type_grouped, tumor_level2, count = n),
                 ratioByGroup = FALSE,
                 showRatioThreshold = 0,
                 labelpositionThreshold = 0.3,
                 showDonutName = TRUE,
                 palette = paletteer::paletteer_d("ggthemes::Classic_Cyclic"))

# Figure 1c -----
meta %>% 
  group_by(metastasis, one_year_progression_free_survival) %>%
  tally() %>%
  na.omit() %>%
  PieDonutCustom(aes(metastasis, one_year_progression_free_survival, count = n),
                 ratioByGroup = TRUE,
                 showRatioThreshold = 0,
                 labelpositionThreshold = 0.3,
                 showDonutName = TRUE,
                 palette = c("#12A2A8FF", "#C7519CFF"))

# Figure 1d -----
df <- meta %>% 
  mutate(day_first_sample = 0) %>%
  mutate(day_last_follow = as.integer(as.POSIXct(date_last_follow) - as.POSIXct(date_first_sample))) 
df %>% 
  pivot_longer(cols = c(day_first_sample, day_last_follow), names_to = "day_type", values_to = "day") %>% 
  mutate(neutropenic_fever = case_when(neutropenic_fever == "Yes" ~ "Neutropenic fever",
                                       neutropenic_fever == "No" ~ "No neutropenic fever",
                                       TRUE ~ neutropenic_fever)) %>%
  mutate(one_year_progression_free_survival = case_when(one_year_progression_free_survival == "Yes" ~ "1yr progression free survival",
                                                        one_year_progression_free_survival == "No" ~ "1yr tumor progression",
                                                        TRUE ~ one_year_progression_free_survival)) %>%
  ggplot(aes(y = factor(study_id, levels = df$study_id[order(df$day_last_follow)]))) +
  geom_line(aes(x = day, color = tumor_type_grouped), size = 2, alpha = 0.3) +
  scale_color_manual(values = c(paletteer::paletteer_d("ggthemes::Classic_Cyclic"))) +
  geom_point(aes(x = day, color = tumor_type_grouped, shape = neutropenic_fever), size = 2) +
  geom_point(aes(x = 365, shape = one_year_progression_free_survival), size = 2, color = "darkgrey", show.legend = FALSE) +
  geom_point(data = subset(df, death == "Yes"), aes(x = day_last_follow), size = 1.5, color = "#B1283AFF", shape = 4, stroke = 2, show.legend = FALSE) +
  scale_shape_manual(values = c("Neutropenic fever" = 16, "No neutropenic fever" = 1, "1yr progression free survival" = 17, "1yr tumor progression" = 2)) +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(x = "Days after baseline samples", y = "Patients", shape = "Outcomes", color = "Tumor type") +
  theme(panel.background = element_blank())
