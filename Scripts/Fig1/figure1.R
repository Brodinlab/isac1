# Script to generate figure1:

# Load packages
library(tidyverse)
library(webr)

# Read functions
source("Scripts/func/PieDonutCustom.R")

# Read metadata table
meta <- read.csv2("Data/isac1_metadata.csv") %>% filter(is.na(subject_id) == FALSE)

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
  labs(y = "Number of patients", x = "Age (Years)") +
  scale_fill_manual(values = c("darkgrey","#A7D3D4FF")) +
  theme(text = element_text(color = "black"), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

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
  mutate(day_last_follow = as.integer(as.POSIXct(date_last_follow) - as.POSIXct(date_first_sample))/365) %>%
  mutate(last_day = day_last_follow) %>% 
  pivot_longer(cols = c(day_first_sample, day_last_follow), names_to = "day_type", values_to = "day")
df$study_id <- factor(df$study_id, levels = unique(df$study_id[order(df$last_day)]))
ggplot(df, aes(y = study_id)) +
  geom_line(aes(x = day, color = tumor_type_grouped), size = 1.8, alpha = 0.3) +
  geom_point(data = subset(df, neutropenic_fever == "Yes"), 
             aes(x = day, color = tumor_type_grouped, fill = tumor_type_grouped), 
             shape = 21, size = 1.5, stroke = 0.2, show.legend = FALSE) +
  geom_point(data = subset(df, neutropenic_fever == "No"), 
             aes(x = day, color = tumor_type_grouped), 
             fill = "white", shape = 21, size = 1.5, stroke = 0.2, show.legend = FALSE) +
  geom_point(data = subset(df, one_year_progression_free_survival == "Yes"), 
             aes(x = 1), 
             color = "black", fill = "black", shape = 24, size = 1.5, stroke = 0.2, show.legend = FALSE) +
  geom_point(data = subset(df, one_year_progression_free_survival == "No"), 
             aes(x = 1), 
             color = "black", fill = NA, shape = 24, size = 1.5, stroke = 0.2, show.legend = FALSE) +
  geom_point(data = subset(df, death == "Yes"), 
             aes(x = last_day), 
             color = "black", shape = 4, size = 1.5, stroke = 0.2, show.legend = FALSE) +
  scale_x_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 1)) +
  scale_color_manual(values = c("#ff9287","#b24502","#d163e6","#00a76c","#006e00","#008cf9","#b80058","#ebac23","#00bbad","#5954d6")) +
  scale_fill_manual(values = c("#ff9287","#b24502","#d163e6","#00a76c","#006e00","#008cf9","#b80058","#ebac23","#00bbad","#5954d6")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(x = "Years after baseline samples", y = "Patients", color = "Tumor type") +
  theme(text = element_text(color = "black"), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
