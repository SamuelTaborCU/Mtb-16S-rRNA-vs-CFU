#######################################
##
## Script name: RMM T95 modeling + association with CFU and 16S rRNA
##
## Purpose of script: (1) Model cure-rate vs treatment duration across RMM studies (Bayesian sigmoidal Emax)
##                    (2) Extract T50/T95 with CIs and save tables
##                    (3) Merge T95 with EOT CFU/16S summaries and generate manuscript plots
##                    (4) Meta-regression + predicted vs observed T95 (and correlation summaries)
##
## Author: Samuel Tabor
##
#######################################

###############################################
# STARTUP

options(stringsAsFactors = FALSE)

library(tidyverse)
library(readxl)
library(rstanemax)
library(broom)
library(grid)     # unit()
library(metafor)

file_path <- "../DataRaw/Manuscript Datasets.xlsx"

plot_dir  <- "../DataProcessed/Plots/5- RMM"
supp_dir  <- file.path(plot_dir, "Supplemental")
stats_dir <- "../DataProcessed/Analysis Data/5- RMM"

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)

###############################################
# IMPORT DATA

RMM <- read_excel(file_path, sheet = "3 RMM studies dataset") %>%
  mutate(
    CFU = as.numeric(CFU),
    Treatment.days = as.numeric(Treatment.days)
  )

###############################################
# DATA PROCESSING (EOT vs Relapse)

MetaDIMS_EOT <- RMM %>%
  filter(Study.stage != "Relapse Assessment") %>%
  filter(!is.na(Adjusted.16S.Original)) %>%
  mutate(Adjusted.16S.Original = round(Adjusted.16S.Original, 0))

MetaDIMS_RLPS <- RMM %>%
  filter(Study.stage == "Relapse Assessment") %>%
  filter(trimws(Culture.status) != "")

###############################################
# RELAPSE SUMMARY (cure_rate)

sum_RLPS_wide <- MetaDIMS_RLPS %>%
  count(Study.name, Group, Treatment.days, Culture.status, name = "N") %>%
  pivot_wider(names_from = Culture.status, values_from = N) %>%
  mutate(
    Negative = ifelse(is.na(Negative), 0, Negative),
    Positive = ifelse(is.na(Positive), 0, Positive),
    cure_rate = Negative / (Positive + Negative),
    Treatment.weeks = Treatment.days / 7
  ) %>%
  arrange(Group, Treatment.days) %>%
  mutate(
    Group = case_when(
      Group == "BPaMZ" & Study.name == "Gates-BC-Rel-01"    ~ "BPaMZ_1",
      Group == "BPaMZ" & Study.name == "TBDA-BC-Rel-01"     ~ "BPaMZ_2",
      Group == "HRZE"  & Study.name == "Gates-BC-Rel-01"    ~ "HRZE_1",
      Group == "HRZE"  & Study.name == "TBDA-BC-Rel-01"     ~ "HRZE_2",
      Group == "HRZE"  & Study.name == "Crush-TB-Substudy"  ~ "HRZE_3",
      Group == "BDOS"  & Study.name == "TBDA-BC-Rel-01"     ~ "BDOS_1",
      Group == "BPaOS" & Study.name == "TBDA-BC-Rel-01"     ~ "BPaOS_1",
      TRUE ~ Group
    )
  ) %>%
  select(Study.name, Group, Treatment.days, Negative, Positive, cure_rate, Treatment.weeks)

###############################################
# MODEL CURE RATES + SAVE PLOT

set.seed(123)

stan_all_RLPS <- stan_emax(
  cure_rate ~ Treatment.weeks,
  data = sum_RLPS_wide,
  e0.fix = NULL,
  emax.fix = NULL,
  gamma.fix = NULL,
  param.cov = list(ec50 = "Group"),
  prior = list(
    e0 = c(0, 1e-3),
    emax = c(1, 1e-3)
  ),
  chains = 2, iter = 5000, seed = 12345
)

plot(stan_all_RLPS, show.ci = FALSE) +
  xlab("Treatment length (weeks)") +
  ylab("% Cured") +
  geom_line(size = 1.5) +
  scale_x_continuous(breaks = c(2,4,6,8,10,12,14,16,18,20)) +
  scale_y_continuous(
    breaks = c(0.0,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00),
    labels = c(0,10,20,30,40,50,60,70,80,90,100)
  ) +
  scale_colour_manual(
    name = "Regimen",
    labels = c(
      "EC50:HRZE_1" = "HRZE",
      "EC50:BPaMZ_1" = "BPaMZ",
      "EC50:PaMZ" = "PaMZ",
      "EC50:BPaL" = "BPaL",
      "EC50:PZM" = "PZM",
      "EC50:BZM" = "BZM",
      "EC50:BZMRb" = "BZMRb",
      "EC50:BDOS_1" = "BDOS",
      "EC50:BPaOS_1" = "BPaOS"
    ),
    values = c(
      "EC50:HRZE_1" = "#E4B811", "EC50:HRZE_2" = "#E4B811", "EC50:HRZE_3" = "#E4B811",
      "EC50:BPaMZ_1" = "red", "EC50:BPaMZ_2" = "red",
      "EC50:PaMZ" = "#00A9FF",
      "EC50:BPaL" = "aquamarine4",
      "EC50:PZM" = "black",
      "EC50:BZM" = "gold4",
      "EC50:BZMRb" = "blue",
      "EC50:BDOS_1" = "forestgreen",
      "EC50:BPaOS_1" = "#C77CFF"
    ),
    breaks = c("EC50:BDOS_1","EC50:BPaL","EC50:BPaMZ_1","EC50:BPaOS_1","EC50:BZM",
               "EC50:BZMRb","EC50:HRZE_1","EC50:PaMZ","EC50:PZM")
  ) +
  scale_fill_manual(
    name = "Regimen",
    labels = c(
      "EC50:HRZE_1" = "HRZE",
      "EC50:BPaMZ_1" = "BPaMZ",
      "EC50:PaMZ" = "PaMZ",
      "EC50:BPaL" = "BPaL",
      "EC50:PZM" = "PZM",
      "EC50:BZM" = "BZM",
      "EC50:BZMRb" = "BZMRb",
      "EC50:BDOS_1" = "BDOS",
      "EC50:BPaOS_1" = "BPaOS"
    ),
    values = c(
      "EC50:HRZE_1" = "#E4B811", "EC50:HRZE_2" = "#E4B811", "EC50:HRZE_3" = "#E4B811",
      "EC50:BPaMZ_1" = "red", "EC50:BPaMZ_2" = "red",
      "EC50:PaMZ" = "#00A9FF",
      "EC50:BPaL" = "aquamarine4",
      "EC50:PZM" = "black",
      "EC50:BZM" = "gold4",
      "EC50:BZMRb" = "blue",
      "EC50:BDOS_1" = "forestgreen",
      "EC50:BPaOS_1" = "#C77CFF"
    ),
    breaks = c("EC50:BDOS_1","EC50:BPaL","EC50:BPaMZ_1","EC50:BPaOS_1","EC50:BZM",
               "EC50:BZMRb","EC50:HRZE_1","EC50:PaMZ","EC50:PZM")
  ) +
  geom_hline(yintercept = 0.5, lty = 3) +
  geom_hline(yintercept = 0.95, lty = 6) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.key.height = unit(1.5, "cm"),
    legend.text = element_text(size = 24, face = "bold"),
    axis.title = element_text(size = 28, color = "black", face = "bold"),
    axis.text.x = element_text(size = 28, color = "black", face = "bold"),
    axis.text.y = element_text(size = 28, color = "black", face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 18)
  )

ggsave(
  filename = file.path(plot_dir, "T95 of the regimens.tiff"),
  width = 9,
  height = 6.75,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "T95 of the regimens.eps"),
  device = cairo_ps,
  width = 9,
  height = 6.75,
  units = "in",
  bg = "white",
  fallback_resolution = 600
)

###############################################
# EXTRACT T50/T95 + SAVE TABLE

sigmoidal_emax_df <- as.data.frame(stan_all_RLPS$stanfit)

ec50_cols <- grep("^ec50\\[", names(sigmoidal_emax_df), value = TRUE)
summary_sigmoidal_emax <- sigmoidal_emax_df[, ec50_cols]

# keep your original naming/order (11 regimens)
names(summary_sigmoidal_emax) <- c(
  "BDOSTBDA","BPaL","BPaMZGates","BPaMZTBDA","BPaOSTBDA",
  "BZM","BZMRb","HRZEGates","HRZECrushTB","PaMZ","PZM"
)

gamma <- quantile(sigmoidal_emax_df$gamma, 0.5)[[1]]

summary_sigmoidal_emax2 <- summary_sigmoidal_emax %>%
  summarise(across(everything(),
                   list(T50 = ~quantile(.x, 0.5), sd = ~sd(.x))))

t50_long <- summary_sigmoidal_emax2[, grepl("T50$", names(summary_sigmoidal_emax2))] %>%
  pivot_longer(everything(), names_to = "Group", values_to = "T50") %>%
  mutate(Group = sub("_T50$", "", Group))

sd_long <- summary_sigmoidal_emax2[, grepl("sd$", names(summary_sigmoidal_emax2))] %>%
  pivot_longer(everything(), names_to = "Group", values_to = "SD") %>%
  mutate(Group = sub("_sd$", "", Group))

merged_df <- left_join(t50_long, sd_long, by = "Group") %>%
  select(Group, SD, T50) %>%
  mutate(
    Group = case_when(
      str_detect(Group, "BPaMZGates")    ~ "BPaMZ-Gates",
      str_detect(Group, "BPaMZTBDA")     ~ "BPaMZ-TBDA",
      str_detect(Group, "HRZEGates")     ~ "HRZE-Gates",
      str_detect(Group, "HRZECrushTB")   ~ "HRZE-Crush-TB",
      str_detect(Group, "BPaOSTBDA")     ~ "BPaOS-TBDA",
      str_detect(Group, "BDOSTBDA")      ~ "BDOS-TBDA",
      TRUE ~ Group
    ),
    Group = factor(
      Group,
      levels = c("HRZE-Gates","HRZE-Crush-TB","BPaMZ-Gates","BPaMZ-TBDA","PaMZ","BPaL",
                 "PZM","BZM","BZMRb","BDOS-TBDA","BPaOS-TBDA")
    )
  ) %>%
  arrange(Group) %>%
  mutate(
    T50_ci_low = T50 - (1.96 * SD),
    T50_ci_up  = T50 + (1.96 * SD),
    T95        = T50 * (95/5)^(1/gamma),
    T95_ci_low = T95 - (1.96 * SD),
    T95_ci_up  = T95 + (1.96 * SD)
  ) %>%
  select(Group, T50, T50_ci_low, T50_ci_up, T95, T95_ci_low, T95_ci_up) %>%
  mutate(
    Study.name = case_when(
      Group %in% c("BPaMZ-Gates","HRZE-Gates","PaMZ","BPaL")          ~ "Gates-BC-Rel-01",
      Group %in% c("BDOS-TBDA","BPaMZ-TBDA","BPaOS-TBDA")             ~ "TBDA-BC-Rel-01",
      Group %in% c("HRZE-Crush-TB","PZM","BZM","BZMRb")               ~ "Crush-TB-Substudy",
      TRUE ~ as.character(Group)
    ),
    Group = as.character(Group),
    Group = case_when(
      Group %in% c("BPaMZ-TBDA","BPaMZ-Gates") ~ "BPaMZ",
      Group %in% c("HRZE-Gates","HRZE-Crush-TB") ~ "HRZE",
      Group == "BDOS-TBDA" ~ "BDOS",
      Group == "BPaOS-TBDA" ~ "BPaOS",
      TRUE ~ Group
    ),
    T95_SE = (T95 - T95_ci_low) / 1.96
  )

write.csv(
  merged_df,
  file.path(stats_dir, "T95,T50 value with CI.csv"),
  row.names = FALSE
)

###############################################
# MERGE T95 WITH EOT MEANS + SAVE DATASETS

sum_EOT <- MetaDIMS_EOT %>%
  group_by(Group, Treatment.days, Study.name) %>%
  summarise(
    mean_CFU = mean(log10(CFU), na.rm = TRUE),
    mean_16S_rRNA = mean(log10(Adjusted.16S.Original), na.rm = TRUE),
    .groups = "drop"
  )

datasetfor_MAX <- merge(MetaDIMS_EOT, merged_df, all = TRUE)

modeling_data <- merge(sum_EOT, merged_df, all = TRUE) %>%
  arrange(Group) %>%
  filter(!Group %in% c("UNTX"))

#write.csv(modeling_data, file.path(stats_dir, "Modeling_data_mean.csv"), row.names = FALSE)
#write.csv(datasetfor_MAX, file.path(stats_dir, "Modeling data.csv"), row.names = FALSE)

# Join PreRx means per study (keeps your behavior)
dfPreRx <- modeling_data %>%
  filter(Group == "PreRx") %>%
  select(Study.name, mean_CFU, mean_16S_rRNA) %>%
  rename(
    mean_CFU_PreRX = mean_CFU,
    mean_16S_rRNA_PreRX = mean_16S_rRNA
  )

modeling_data <- modeling_data %>%
  left_join(dfPreRx, by = "Study.name")

New_modeling_data <- modeling_data %>%
  filter(Treatment.days %in% c(7,14,21,28,42,56)) %>%
  mutate(Treatment.days = as.numeric(Treatment.days)) %>%
  arrange(Treatment.days)

modeling_data_28 <- modeling_data %>% filter(Treatment.days == 28)

###############################################
# PLOTS (unchanged output logic/files)

CFU <- ggplot(data = New_modeling_data, aes(T95, mean_CFU, col=Group))+
  geom_point(size=4)+
  scale_colour_manual(values = c("PreRx"="grey50","HRZE"= "#E4B811","PaMZ"= "#00A9FF","BDOS"="forestgreen",
                                 "BPaOS"= "#C77CFF","BPaL"="aquamarine4","PZM"="black",
                                 "BZM"="gold4","BZMRb"="blue","BPaMZ"="red"),) +
  scale_y_continuous(limits = c(0,7), breaks=c(0,1,2,3,4,5,6,7))+
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  xlab("T95 In Weeks")+
  ylab("CFU(log 10)")+ facet_wrap(~Treatment.days) +
  ggtitle("CFU vs T95 in different timepoints")+
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.key.height = unit(1.15, "cm"),
        legend.text = element_text(size = 20, face = "bold")) +
  theme(axis.title = element_text(size = 28, color = "black", face = "bold")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(axis.text.x = element_text(size = 28, color = "black", face = "bold"),
        axis.text.y = element_text(size = 28, color = "black", face = "bold")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5, size = 18))

ggsave(
  filename = file.path(supp_dir, "CFU VS T95 at different timepoints.tiff"),
  plot = CFU,
  width = 9,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = file.path(supp_dir, "CFU VS T95 at different timepoints.eps"),
  plot = CFU,
  device = cairo_ps,
  width = 9,
  height = 6,
  units = "in",
  bg = "white",
  fallback_resolution = 600
)

CFU28 <- ggplot(data = modeling_data_28 , aes(T95, mean_CFU, col=Group))+
  geom_point(size=4)+
  scale_colour_manual(values = c("PreRx"="grey50","HRZE"= "#E4B811","PaMZ"= "#00A9FF","BDOS"="forestgreen",
                                 "BPaOS"= "#C77CFF","BPaL"="aquamarine4","PZM"="black",
                                 "BZM"="gold4","BZMRb"="blue","BPaMZ"="red"),) +
  scale_y_continuous(limits = c(0,6.5), breaks=c(0,1,2,3,4,5,6,7))+
  xlab("T95 In Weeks")+
  ylab("CFU (Log10)")+
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.key.height = unit(1.15, "cm"),
        legend.text = element_text(size = 20, face = "bold")) +
  theme(axis.title = element_text(size = 28, color = "black", face = "bold")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(axis.text.x = element_text(size = 28, color = "black", face = "bold"),
        axis.text.y = element_text(size = 28, color = "black", face = "bold")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5, size = 18))

ggsave(
  filename = file.path(supp_dir, "CFU VS T95 day 28.tiff"),
  plot = CFU28,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = file.path(supp_dir, "CFU VS T95 day 28.eps"),
  plot = CFU28,
  device = cairo_ps,
  width = 7,
  height = 6,
  units = "in",
  bg = "white",
  fallback_resolution = 600
)

A16S.rRNA <- ggplot(data = New_modeling_data, aes( T95, mean_16S_rRNA, col=Group))+
  geom_point(size=4)+
  scale_colour_manual(values = c("PreRx"="grey50","HRZE"= "#E4B811","PaMZ"= "#00A9FF","BDOS"="forestgreen",
                                 "BPaOS"= "#C77CFF","BPaL"="aquamarine4","PZM"="black",
                                 "BZM"="gold4","BZMRb"="blue","BPaMZ"="red"),) +
  scale_y_continuous(limits = c(4,9.5), breaks=c(0,1,2,3,4,5,6,7,8,9,10))+
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  xlab("T95 In weeks")+
  ylab("16S.rRNA(log10)")+ facet_wrap(~Treatment.days) +
  ggtitle("16S rRNA vs T95 in different timepoints")+
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.key.height = unit(1.15, "cm"),
        legend.text = element_text(size = 20, face = "bold")) +
  theme(axis.title = element_text(size = 28, color = "black", face = "bold")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(axis.text.x = element_text(size = 28, color = "black", face = "bold"),
        axis.text.y = element_text(size = 28, color = "black", face = "bold")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5, size = 18))

ggsave(
  filename = file.path(supp_dir, "16S rRNA VS T95 at different timepoints.tiff"),
  plot = A16S.rRNA,
  width = 9,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = file.path(supp_dir, "16S rRNA VS T95 at different timepoints.eps"),
  plot = A16S.rRNA,
  device = cairo_ps,
  width = 9,
  height = 6,
  units = "in",
  bg = "white",
  fallback_resolution = 600
)
A16S.rRNA28 <- ggplot(data = modeling_data_28, aes( T95, mean_16S_rRNA, col=Group))+
  geom_point(size=4)+
  scale_colour_manual(values = c("PreRx"="grey50","HRZE"= "#E4B811","PaMZ"= "#00A9FF","BDOS"="forestgreen",
                                 "BPaOS"= "#C77CFF","BPaL"="aquamarine4","PZM"="black",
                                 "BZM"="gold4","BZMRb"="blue","BPaMZ"="red"),) +
  scale_y_continuous(limits = c(4,9.5), breaks=c(0,1,2,3,4,5,6,7,8,9,10))+
  xlab("T95 In weeks")+
  ylab("16S rRNA (Log10)")+
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.key.height = unit(1.15, "cm"),
        legend.text = element_text(size = 20, face = "bold")) +
  theme(axis.title = element_text(size = 28, color = "black", face = "bold")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(axis.text.x = element_text(size = 28, color = "black", face = "bold"),
        axis.text.y = element_text(size = 28, color = "black", face = "bold")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5, size = 18))

ggsave(
  filename = file.path(supp_dir, "16S rRNA VS T95 day 28.tiff"),
  plot = A16S.rRNA28,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = file.path(supp_dir, "16S rRNA VS T95 day 28.eps"),
  plot = A16S.rRNA28,
  device = cairo_ps,
  width = 7,
  height = 6,
  units = "in",
  bg = "white",
  fallback_resolution = 600
)
###############################################
# META-REGRESSION 

dataRaw_path <- "../DataRaw/Modeling data.csv"
dataRaw <- read.csv(dataRaw_path)

datasetfor_MAX <- dataRaw

t95Values <- datasetfor_MAX |>
  filter(Treatment.days == 28, Study.stage == "On Treatment") |>
  select(STUDY = Study.name, GROUP = Group, TREATMENT_DAYS = Treatment.days, T95,
         T95_CI_LOW = T95_ci_low, T95_CI_UP = T95_ci_up) |>
  distinct() |>
  select(-TREATMENT_DAYS)

dataClean <- datasetfor_MAX |>
  group_by(Study.name, Group, Treatment.days) |>
  summarize(
    MEAN_A16S = mean(Adjusted.16S.Original, na.rm = TRUE),
    MEAN_LOG_A16S = mean(log10(Adjusted.16S.Original), na.rm = TRUE),
    MEAN_CFU = mean(CFU, na.rm = TRUE),
    MEAN_LOG_CFU = mean(log10(CFU), na.rm = TRUE),
    N = n(),
    .groups = "drop"
  ) |>
  rename(STUDY = Study.name, GROUP = Group, TREATMENT_DAYS = Treatment.days)

dataCleanDay0 <- dataClean |>
  filter(TREATMENT_DAYS == 0) |>
  select(
    STUDY,
    MEAN_A16S_BL = MEAN_A16S,
    MEAN_LOG_A16S_BL = MEAN_LOG_A16S,
    MEAN_CFU_BL = MEAN_CFU,
    MEAN_LOG_CFU_BL = MEAN_LOG_CFU,
    N_BL = N
  )

dataClean <- dataClean |>
  filter(TREATMENT_DAYS == 28) |>
  select(
    STUDY, GROUP,
    MEAN_A16S_D28 = MEAN_A16S,
    MEAN_CFU_D28 = MEAN_CFU,
    MEAN_LOG_A16S_D28 = MEAN_LOG_A16S,
    MEAN_LOG_CFU_D28 = MEAN_LOG_CFU,
    N_D28 = N
  ) |>
  left_join(dataCleanDay0, by = "STUDY") |>
  left_join(t95Values, by = c("STUDY", "GROUP")) |>
  mutate(
    MEAN_A16S_DIFF = MEAN_A16S_D28 - MEAN_A16S_BL,
    MEAN_LOG_A16S_DIFF = MEAN_LOG_A16S_D28 - MEAN_LOG_A16S_BL,
    MEAN_CFU_DIFF = MEAN_CFU_D28 - MEAN_CFU_BL,
    MEAN_LOG_CFU_DIFF = MEAN_LOG_CFU_D28 - MEAN_LOG_CFU_BL,
    T95_SE = (T95 - T95_CI_LOW) / 1.96
  )

blAnd28FitCfu <- metafor::rma(
  yi = T95, sei = T95_SE,
  mods = ~ MEAN_LOG_CFU_BL + MEAN_LOG_CFU_DIFF,
  data = dataClean, method = "ML"
)

blAnd28Fit16s <- metafor::rma(
  yi = T95, sei = T95_SE,
  mods = ~ MEAN_LOG_A16S_BL + MEAN_LOG_A16S_DIFF,
  data = dataClean, method = "ML"
)

blAnd28FitBoth <- metafor::rma(
  yi = T95, sei = T95_SE,
  mods = ~ MEAN_LOG_A16S_BL + MEAN_LOG_A16S_DIFF + MEAN_LOG_CFU_BL + MEAN_LOG_CFU_DIFF,
  data = dataClean, method = "ML"
)

nullFit <- metafor::rma(yi = T95, sei = T95_SE, mods = ~ 1, data = dataClean, method = "ML")

a16sAnova <- metafor::anova.rma(blAnd28FitBoth, blAnd28Fit16s)
cfuAnova  <- metafor::anova.rma(blAnd28FitBoth, blAnd28FitCfu)

a16sNullAnova <- metafor::anova.rma(nullFit, blAnd28Fit16s)
cfuNullAnova  <- metafor::anova.rma(nullFit, blAnd28FitCfu)
bothNullAnova <- metafor::anova.rma(nullFit, blAnd28FitBoth)

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

hypoThreeTable <- data.frame(
  Variables = c(
    "BL CFU + 28-day change CFU",
    "BL 16S rRNA + 28-day change 16S rRNA",
    "BL CFU + 28-day change CFU + BL 16S rRNA + 28-day change 16S rRNA"
  ),
  R2 = c(blAnd28FitCfu$R2, blAnd28Fit16s$R2, blAnd28FitBoth$R2),
  AIC = c(AIC(blAnd28FitCfu), AIC(blAnd28Fit16s), AIC(blAnd28FitBoth)),
  RMSE = c(
    sqrt(mean(residuals(blAnd28FitCfu)^2)),
    sqrt(mean(residuals(blAnd28Fit16s)^2)),
    sqrt(mean(residuals(blAnd28FitBoth)^2))
  ),
  lrtPVal_raw = c(cfuAnova$pval, a16sAnova$pval, NA),
  lrtNullPval_raw = c(cfuNullAnova$pval, a16sNullAnova$pval, bothNullAnova$pval),
  check.names = FALSE
)

hypoThreeTable[["Pseudo R^2"]] <- round(hypoThreeTable$R2 / 100, 2)
hypoThreeTable$AIC <- round(hypoThreeTable$AIC, 2)
hypoThreeTable$RMSE <- round(hypoThreeTable$RMSE, 2)

hypoThreeTable[["LRT P-value"]] <- ifelse(is.na(hypoThreeTable$lrtPVal_raw), "Ref.", fmt_p(hypoThreeTable$lrtPVal_raw))
hypoThreeTable[["LRT P-value (null comparison)"]] <- fmt_p(hypoThreeTable$lrtNullPval_raw)

hypoThreeTable <- hypoThreeTable[, c("Variables","Pseudo R^2","AIC","RMSE","LRT P-value","LRT P-value (null comparison)")]

write.csv(hypoThreeTable, file.path(stats_dir, "Output from meta regression modeling of T95.csv"), row.names = FALSE)

predictedDf <- data.frame(
  T95         = dataClean$T95,
  Group       = dataClean$GROUP,
  pred_model1 = predict(blAnd28FitCfu)$pred,
  pred_model2 = predict(blAnd28Fit16s)$pred,
  pred_model3 = predict(blAnd28FitBoth)$pred
)

get_correlation_results <- function(actual, predicted) {
  r <- cor.test(actual, predicted)
  data.frame(
    Correlation_Coefficient = unname(r$estimate),
    CI_Lower = r$conf.int[1],
    CI_Upper = r$conf.int[2],
    P_Value = r$p.value
  )
}

correlation_summary <- rbind(
  model1 = get_correlation_results(predictedDf$T95, predictedDf$pred_model1),
  model2 = get_correlation_results(predictedDf$T95, predictedDf$pred_model2),
  model3 = get_correlation_results(predictedDf$T95, predictedDf$pred_model3)
)

write.csv(correlation_summary, file.path(stats_dir, "Correlation actual T95 vs Pridicted T95.csv"), row.names = FALSE)

cols_map <- c(
  "PreRx"="grey50","HRZE"= "#E4B811","PaMZ"= "#00A9FF","BDOS"="forestgreen",
  "BPaOS"= "#C77CFF","BPaL"="aquamarine4","PZM"="black",
  "BZM"="gold4","BZMRb"="blue","BPaMZ"="red"
)

make_model_plot <- function(df, yvar, ylab_text, limits = c(3, 22)) {
  ggplot(df, aes(x = T95, y = .data[[yvar]], col = Group)) +
    geom_point(size = 4) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 1) +
    scale_colour_manual(values = cols_map) +
    xlab("T95 In Weeks") +
    ylab(ylab_text) +
    scale_x_continuous(limits = limits, breaks = c(0,5,10,15,20,25,30)) +
    scale_y_continuous(limits = limits, breaks = c(0,5,10,15,20,25,30)) +
    coord_equal(xlim = limits, ylim = limits, expand = FALSE) +
    theme_bw() +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      legend.key.height = unit(1.15, "cm"),
      legend.text = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 28, color = "black", face = "bold"),
      axis.text.x = element_text(size = 28, color = "black", face = "bold"),
      axis.text.y = element_text(size = 28, color = "black", face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      plot.title = element_text(hjust = 0.5, size = 18)
    )
}

model1 <- make_model_plot(predictedDf, "pred_model1", "Model 1 T95 Prediction")
model2 <- make_model_plot(predictedDf, "pred_model2", "Model 2 T95 Prediction")
model3 <- make_model_plot(predictedDf, "pred_model3", "Model 3 T95 Prediction")

ggsave(
  filename = file.path(plot_dir, "Model 1 T95 predicted vs observed.tiff"),
  plot = model1,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "Model 1 T95 predicted vs observed.eps"),
  plot = model1,
  device = cairo_ps,
  width = 7,
  height = 6,
  units = "in",
  bg = "white",
  fallback_resolution = 600
)

ggsave(
  filename = file.path(plot_dir, "Model 2 T95 predicted vs observed.tiff"),
  plot = model2,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "Model 2 T95 predicted vs observed.eps"),
  plot = model2,
  device = cairo_ps,
  width = 7,
  height = 6,
  units = "in",
  bg = "white",
  fallback_resolution = 600
)

ggsave(
  filename = file.path(plot_dir, "Model 3 T95 predicted vs observed.tiff"),
  plot = model3,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "Model 3 T95 predicted vs observed.eps"),
  plot = model3,
  device = cairo_ps,
  width = 7,
  height = 6,
  units = "in",
  bg = "white",
  fallback_resolution = 600
)
#############################################
############## The End ######################
#############################################