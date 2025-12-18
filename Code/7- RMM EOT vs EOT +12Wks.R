#######################################
##
## Script name: RMM EOT vs 12-week relapse assessment summary (Table 4)
##
## Purpose of script: (1) Subset RMM dataset to EOT (On Treatment) and relapse assessment (12 weeks post-treatment)
##                    (2) Compute n and % culture-negative at EOT and relapse assessment
##                    (3) Compute median log10 CFU and median log10 16S rRNA at EOT and relapse assessment
##                    (4) Calculate change from EOT to relapse assessment and save a clean CSV table
##
## Author: Samuel Tabor
##
#######################################

###############################################
####################STARTUP####################
###############################################

# Import libraries
# Set options
# Define variables

options(stringsAsFactors = FALSE)

library(readxl)
library(dplyr)
library(tidyr)

file_path <- "../DataRaw/Manuscript Datasets.xlsx"
out_csv   <- "../DataProcessed/Analysis Data/4- Rebound/RMM_Table4_ordered.csv"

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)

###############################################
##################IMPORT DATA##################
###############################################

rmm <- read_excel(file_path, sheet = "3 RMM studies dataset")

###############################################
#################DATA PROCESSING###############
###############################################

rmm_clean <- rmm %>%
  filter(!is.na(CFU), !is.na(Adjusted.16S.Original)) %>%
  mutate(
    CFU = as.numeric(CFU),
    Adjusted.16S.Original = as.numeric(Adjusted.16S.Original),
    log10CFU  = if_else(CFU > 0, log10(CFU), NA_real_),
    log10_16S = if_else(Adjusted.16S.Original > 0, log10(Adjusted.16S.Original), NA_real_)
  )

eot <- rmm_clean %>% filter(Study.stage == "On Treatment")
ra  <- rmm_clean %>% filter(Study.stage == "Relapse Assessment")

eot_counts <- eot %>%
  group_by(Group, Treatment.days) %>%
  summarise(
    n_EOT        = n(),
    n_neg_EOT    = sum(Culture.status == "Negative", na.rm = TRUE),
    frac_neg_EOT = n_neg_EOT / n_EOT,
    .groups      = "drop"
  )

ra_counts <- ra %>%
  group_by(Group, Treatment.days) %>%
  summarise(
    n_RA        = n(),
    n_neg_RA    = sum(Culture.status == "Negative", na.rm = TRUE),
    frac_neg_RA = n_neg_RA / n_RA,
    .groups     = "drop"
  )

eot_med <- eot %>%
  group_by(Group, Treatment.days) %>%
  summarise(
    Median_CFU_EOT = median(log10CFU, na.rm = TRUE),
    Median_16S_EOT = median(log10_16S, na.rm = TRUE),
    .groups        = "drop"
  )

ra_med <- ra %>%
  group_by(Group, Treatment.days) %>%
  summarise(
    Median_CFU_RA  = median(log10CFU, na.rm = TRUE),
    Median_16S_RA  = median(log10_16S, na.rm = TRUE),
    .groups        = "drop"
  )

table4 <- eot_med %>%
  inner_join(ra_med, by = c("Group", "Treatment.days")) %>%
  left_join(eot_counts, by = c("Group", "Treatment.days")) %>%
  left_join(ra_counts,  by = c("Group", "Treatment.days")) %>%
  mutate(
    Duration_weeks        = Treatment.days / 7,
    Change_CFU_log10      = Median_CFU_RA  - Median_CFU_EOT,
    Change_16S_log10      = Median_16S_RA  - Median_16S_EOT,
    Percent_negative_EOT  = round(frac_neg_EOT * 100, 1),
    Percent_negative_RA   = round(frac_neg_RA * 100, 1)
  )

table4_formatted <- table4 %>%
  arrange(Group, Duration_weeks) %>%
  transmute(
    Regimen                  = Group,
    Duration_weeks           = as.numeric(Duration_weeks),
    
    n_EOT                    = n_EOT,
    n_negative_EOT           = n_neg_EOT,
    Percent_negative_EOT     = Percent_negative_EOT,
    
    n_RA                     = n_RA,
    n_negative_RA            = n_neg_RA,
    Percent_negative_RA      = Percent_negative_RA,
    
    Median_CFU_EOT           = round(Median_CFU_EOT, 2),
    Median_CFU_EOT_plus12wks = round(Median_CFU_RA, 2),
    Change_CFU_log10         = round(Change_CFU_log10, 2),
    
    Median_16S_EOT           = round(Median_16S_EOT, 2),
    Median_16S_EOT_plus12wks = round(Median_16S_RA, 2),
    Change_16S_log10         = round(Change_16S_log10, 2)
  )

table4_filtered <- table4_formatted %>%
  filter((100 - Percent_negative_EOT) > 33)

regimen_order <- c("BPaL", "BPaMZ", "HRZE", "BDOS", "BPaOS")

table4_ordered <- table4_filtered %>%
  mutate(Regimen = factor(Regimen, levels = regimen_order)) %>%
  arrange(Regimen, Duration_weeks)

table4_ordered
write.csv(table4_ordered, out_csv, row.names = FALSE)
