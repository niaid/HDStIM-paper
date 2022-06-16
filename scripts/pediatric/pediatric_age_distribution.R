#!/usr/bin/env Rscript --vanilla

# PURPOSE: Generate stats and age distribution plots for the pediatric dataset. 

library(tidyverse)

results_folder <- file.path("results", "pediatric")
figures_folder <- file.path("figures", "pediatric")

# Selected data. 
selected_data <- readRDS(file.path(results_folder, "selected_data_all.rds"))
state_markers <- selected_data$state_markers

# Panel data.
panel <- read_tsv(file.path("meta", "pediatric_stim_panel.txt"))

# Args file.
args <- read_tsv(file.path("meta", "pediatric_fcs_info.txt"))

# Number of case and control subjects.
case_subj <- dplyr::filter(args, condition == "Case")
no_case <- case_subj %>% dplyr::select(patient_id) %>% unique

control_subj <- dplyr::filter(args, condition == "Control")
no_control <- control_subj %>% dplyr::select(patient_id) %>% unique

# Age range and gender.
case_range <- selected_data$response_mapping_main %>% 
  dplyr::filter(condition == "Case") %>%
  dplyr::select(age_at_sample) %>% unique

min(case_range)
max(case_range)

control_range <- selected_data$response_mapping_main %>% 
  dplyr::filter(condition == "Control") %>%
  dplyr::select(age_at_sample) %>% unique

min(control_range)
max(control_range)

case_gender <- selected_data$response_mapping_main %>% 
  dplyr::filter(condition == "Case") %>%
  dplyr::select(patient_id, gender) %>%
  unique()
table(case_gender$gender)

control_gender <- selected_data$response_mapping_main %>% 
  dplyr::filter(condition == "Control") %>%
  dplyr::select(patient_id, gender) %>%
  unique()
table(control_gender$gender)

all_male <- selected_data$response_mapping_main %>% 
  dplyr::filter(gender == "M") %>%
  dplyr::select(patient_id, gender, age_at_sample) %>%
  unique()
table(all_male$gender)
min(all_male$age_at_sample)
max(all_male$age_at_sample)

all_female <- selected_data$response_mapping_main %>% 
  dplyr::filter(gender == "F") %>%
  dplyr::select(patient_id, gender, age_at_sample) %>%
  unique()
table(all_female$gender)
min(all_female$age_at_sample)
max(all_female$age_at_sample)

# Plot a histogram. 
all_subj <- rbind(all_male, all_female)
hp <- ggplot(all_subj, aes(x = age_at_sample, fill = gender)) +
  geom_histogram(binwidth = 1, color = "black", alpha=0.9) +
  labs(x = "Age", y = "Count", fill = "Gender") +
  scale_x_continuous(name = "Age", breaks = c(2:16)) +
  theme(axis.text = element_text(size = 16)) +
  theme(axis.title = element_text(size = 18)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.title = element_text(size = 14))

ggsave("pediatric_age_distribution.png", plot = hp, path = figures_folder, width = 6, height = 4, dpi = 300)
