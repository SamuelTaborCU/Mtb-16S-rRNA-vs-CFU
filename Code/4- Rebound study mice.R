#######################################
##
## Script name: Rebound study
##
## Purpose of script: Generate manuscript plots and summaries for the rebound study:
##                    (A) longitudinal CFU and 16S rRNA time courses (Weeks 2 and 4 HRZE)
##                    (B) CFU vs 16S rRNA reductions from baseline by treatment duration (WK4 and WK8)
##                    (C) statistical summaries
##
## Author: Samuel Tabor
##
#######################################

###############################################
####################STARTUP####################
###############################################
# Import libraries / Set options / Define variables

options(stringsAsFactors = FALSE)

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(broom)
library(jpeg)
library(grid)

file_path <- "../DataRaw/Manuscript Datasets.xlsx"
plot_dir  <- "../DataProcessed/Plots/4- Rebound"
stats_dir <- "../DataProcessed/Analysis Data/4- Rebound"

dir.create(plot_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)

###############################################
##################IMPORT DATA##################
###############################################

rna_rebound <- read_excel(file_path, sheet = "Rebound")

###############################################
#################DATA PROCESSING###############
###############################################

# NOTE: Data are assumed clean (no "<" strings; no NEGATIVE-culture adjustment needed)

rna_rebound_1 <- rna_rebound %>%
  mutate(
    Days.since.Tx.start = as.numeric(Days.since.Tx.start),
    Treatment.days      = as.numeric(Treatment.days),
    Washout.period      = as.numeric(Washout.period),
    Washout.period      = ifelse(is.na(Washout.period), 0, Washout.period),
    rebound_date        = Days.since.Tx.start - Treatment.days,
    
    CFU                 = as.numeric(CFU),
    Adjusted.16S.Original = as.numeric(Adjusted.16S.Original),
    
    Group2 = ifelse(Washout.period %in% c(0, 1), "PreRx", as.character(Group))
  )

# N per Group x Day (wide)
sum_rna_rebound_1_wide <- rna_rebound_1 %>%
  group_by(Group, Days.since.Tx.start) %>%
  summarise(N = n(), .groups = "drop") %>%
  pivot_wider(id_cols = Group, names_from = Days.since.Tx.start, values_from = N) %>%
  as.data.frame()

print(sum_rna_rebound_1_wide)

# Keep only groups used in rebound plots
selected_groups <- c("UNTX", "HRZE 2 WK DOSE", "HRZE 4 WK DOSE", "HRZE 8 WK DOSE")

rna_rebound_2 <- rna_rebound_1 %>%
  mutate(
    Group = factor(Group, levels = selected_groups, labels = selected_groups)
  ) %>%
  filter(Group %in% selected_groups) %>%
  mutate(
    Group2 = factor(
      Group2,
      levels = c("PreRx", "HRZE 2 WK DOSE", "HRZE 4 WK DOSE", "HRZE 8 WK DOSE"),
      labels = c("PreRx", "HRZE 2 WK DOSE", "HRZE 4 WK DOSE", "HRZE 8 WK DOSE")
    ),
    
    Group3 = as.character(Group2),
    Group3 = ifelse(Group == "UNTX", "PreRx", Group3),
    Group3 = ifelse(Washout.period == 1, "On Treatment", Group3),
    Group3 = factor(
      Group3,
      levels = c("PreRx", "On Treatment", "HRZE 2 WK DOSE", "HRZE 4 WK DOSE", "HRZE 8 WK DOSE"),
      labels = c("PreRx", "On Treatment", "HRZE 2 WK DOSE", "HRZE 4 WK DOSE", "HRZE 8 WK DOSE")
    )
  )

median_by <- function(df, grp_col, y_col) {
  df %>%
    filter(!is.na(.data[[y_col]])) %>%
    group_by(grp = .data[[grp_col]], Days.since.Tx.start) %>%
    summarise(val = median(.data[[y_col]], na.rm = TRUE), .groups = "drop") %>%
    rename(!!grp_col := grp, !!y_col := val) %>%
    as.data.frame()
}

###############################################
###############PLOTS AND TABLES################
###############################################
# Week 2 HRZE longitudinal plots (Off-treatment = HRZE 2)

rna_w2 <- rna_rebound_2 %>%
  filter(Group3 != "HRZE 4 WK DOSE" & Group3 != "HRZE 8 WK DOSE")

rna_w2_CFU <- rna_w2 %>%
  filter(!is.na(CFU)) %>%
  filter(Group %in% c("UNTX", "HRZE 2 WK DOSE"))

rna_w2_16S <- rna_w2 %>%
  filter(!is.na(Adjusted.16S.Original)) %>%
  filter(Group %in% c("UNTX", "HRZE 2 WK DOSE"))

sum_CFU_w2  <- median_by(rna_w2_CFU, "Group",  "CFU")
sum_CFU2_w2 <- median_by(rna_w2_CFU, "Group2", "CFU")

sum_16S_w2  <- median_by(rna_w2_16S, "Group",  "Adjusted.16S.Original")
sum_16S2_w2 <- median_by(rna_w2_16S, "Group2", "Adjusted.16S.Original") %>%
  dplyr::filter(Days.since.Tx.start %in% c(0, 12))

# ---- CFU plot (Week 2) ----
longitudinal_CFU_w2 <- ggplot(
  rna_w2_CFU,
  aes(
    Days.since.Tx.start, CFU,
    color = factor(Group2, levels = c("PreRx", "HRZE 2 WK DOSE", "HRZE 4 WK DOSE"))
  )
) +
  geom_jitter(width = 0.0, size = 5) +
  geom_line(data = sum_CFU_w2, aes(Days.since.Tx.start, CFU, group = Group, col = Group), size = 2) +
  geom_line(
    data = subset(sum_CFU2_w2, Days.since.Tx.start <= 12),
    aes(Days.since.Tx.start, CFU, group = Group2, col = Group2),
    size = 2
  ) +
  scale_color_manual(values = c("PreRx" = "brown", "On Treatment" = "brown", "HRZE 2 WK DOSE" = "grey25"),
                     labels = c("PreRx" = "On Treatment", "HRZE 2 WK DOSE" = "Off Treatment")) +
  scale_y_log10(limits = c(7e+03, 1e+08), breaks = c(1e+01,1e+02,1e+03,1e+04,1e+05,1e+06,1e+07,1e+08),
                labels = c(1,2,3,4,5,6,7,8)) +
  scale_x_continuous(limits = c(0,56), breaks = c(0,14,28,42,56,70,84,98,112),
                     labels = c("0", "2","4","6","8","10 weeks","12 weeks","14 weeks","16 weeks")) +
  geom_vline(xintercept = 12, linetype = "dashed", color = "red3") +
  geom_hline(yintercept = 726562.5, linetype = "dashed", color = "grey50") +
  annotate("text", x = 5, y = 0, label = "Treatment",
           hjust = 0.5, vjust = -27.5, color = "black", size = 5, fontface = "bold") +
  annotate("text", x = 5, y = 0, label = " Recovery Phase",
           hjust = -1.75, vjust = -27.5, color = "black", size = 5, fontface = "bold") +
  labs(y = "CFU (Log10)", x = "Weeks", colour = "Treatment Status") +
  theme_bw() +
  theme(legend.position = "", axis.title = element_text(family = "", color = "black", size = 18)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(axis.text.x = element_text(size = 28, color = "black", face = "bold"),
        axis.text.y = element_text(size = 28, color = "black", face = "bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5, size = 18))

longitudinal_CFU_w2
ggsave(file.path(plot_dir, "Rebound longitudinal CFU Week 2 HRZE.jpg"),
       longitudinal_CFU_w2, width = 7, height = 5)

# ---- 16S plot (Week 2) ----
longitudinal_16S_w2 <- ggplot(
  rna_w2_16S,
  aes(Days.since.Tx.start, Adjusted.16S.Original, col = Group2)
) +
  geom_jitter(width = 0.0, size = 5, shape = 17) +
  geom_line(data = sum_16S_w2, aes(Days.since.Tx.start, Adjusted.16S.Original, group = Group, col = Group), size = 2) +
  geom_line(data = sum_16S2_w2, aes(Days.since.Tx.start, Adjusted.16S.Original, group = Group2, col = Group2), size = 2) +
  scale_y_log10(limits = c(6e+06, 1e+10), breaks = c(1e+04,1e+05,1e+06,1e+07,1e+8,1e+09,1e+10),
                labels = c(4,5,6,7,8,9,10)) +
  annotate("text", x = 5, y = 0, label = "Treatment",
           hjust = 0.5, vjust = -27.5, color = "black", size = 5, fontface = "bold") +
  annotate("text", x = 5, y = 0, label = " Recovery Phase",
           hjust = -1.75, vjust = -27.5, color = "black", size = 5, fontface = "bold") +
  scale_x_continuous(limits = c(0,56), breaks = c(0,14,28,42,56,70,84,98,112),
                     labels = c("0", "2","4","6","8","10 weeks","12 weeks","14 weeks","16 weeks")) +
  scale_color_manual(values = c("PreRx" = "brown", "On Treatment" = "brown", "HRZE 2 WK DOSE" = "grey25"),
                     labels = c("PreRx" = "On Treatment", "HRZE 2 WK DOSE" = "Off Treatment")) +
  labs(y = "16S rRNA (log10)", x = "Weeks", colour = "Treatment Status") +
  geom_hline(yintercept = 280111899, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 12, linetype = "dashed", color = "red3") +
  theme_bw() +
  theme(legend.position = "", axis.title = element_text(family = "", color = "black", size = 18)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(axis.text.x = element_text(size = 28, color = "black", face = "bold"),
        axis.text.y = element_text(size = 28, color = "black", face = "bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5, size = 18))

longitudinal_16S_w2
ggsave(file.path(plot_dir, "Rebound longitudinal 16SrRNA Week 2 HRZE.jpg"),
       longitudinal_16S_w2, width = 7, height = 5)

plot_grid(longitudinal_CFU_w2, longitudinal_16S_w2, align = "hv", nrow = 1)

###############################################
###############PLOTS AND TABLES################
###############################################
# Week 4 HRZE longitudinal plots (Off-treatment = HRZE 4)

rna_w4 <- rna_rebound_2 %>%
  filter(Group3 != "HRZE 2 WK DOSE" & Group3 != "HRZE 8 WK DOSE")

rna_w4_CFU <- rna_w4 %>% filter(!is.na(CFU))
rna_w4_16S <- rna_w4 %>% filter(!is.na(Adjusted.16S.Original))

sum_CFU_w4  <- median_by(rna_w4_CFU, "Group",  "CFU")
sum_CFU2_w4 <- median_by(rna_w4_CFU, "Group2", "CFU")

sum_16S_w4  <- median_by(rna_w4_16S, "Group",  "Adjusted.16S.Original")
sum_16S2_w4 <- median_by(rna_w4_16S, "Group2", "Adjusted.16S.Original") %>%
  dplyr::filter(Days.since.Tx.start %in% c(0, 12, 26))

# ---- CFU plot (Week 4) ----
longitudinal_CFU_w4 <- ggplot(
  rna_w4_CFU,
  aes(
    Days.since.Tx.start, CFU,
    color = factor(Group2, levels = c("PreRx", "HRZE 2 WK DOSE", "HRZE 4 WK DOSE", "HRZE 8 WK DOSE"))
  )
) +
  geom_jitter(width = 0.0, size = 4) +
  geom_line(data = sum_CFU_w4, aes(Days.since.Tx.start, CFU, group = Group, col = Group), size = 2) +
  geom_line(data = sum_CFU2_w4, aes(Days.since.Tx.start, CFU, group = Group2, col = Group2), size = 2) +
  annotate("text", x = 5, y = 0, label = "Treatment",
           hjust = 0, vjust = -27.5, color = "black", size = 5, fontface = "bold") +
  annotate("text", x = 5, y = 0, label = " Recovery Phase",
           hjust = -2, vjust = -27.5, color = "black", size = 5, fontface = "bold") +
  scale_color_manual(values = c("PreRx" = "brown", "On Treatment" = "brown", "HRZE 4 WK DOSE" = "grey25"),
                     labels = c("PreRx" = "On Treatment", "HRZE 4 WK DOSE" = "Off Treatment")) +
  scale_y_log10(limits = c(7e+03, 1e+08), breaks = c(1e+01,1e+02,1e+03,1e+04,1e+05,1e+06,1e+07,1e+08),
                labels = c(1,2,3,4,5,6,7,8)) +
  scale_x_continuous(limits = c(0,56), breaks = c(0,14,28,42,56,70,84,98,112),
                     labels = c("0", "2","4","6","8","10","12 weeks","14 weeks","16 weeks")) +
  geom_vline(xintercept = 26, linetype = "dashed", color = "red3") +
  geom_hline(yintercept = 257812.50, linetype = "dashed", color = "grey50") +
  labs(y = "CFU (Log10)", x = "Weeks", colour = "Treatment Status") +
  theme_bw() +
  theme(legend.position = "", axis.title = element_text(family = "", color = "black", size = 18)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(axis.text.x = element_text(size = 28, color = "black", face = "bold"),
        axis.text.y = element_text(size = 28, color = "black", face = "bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5, size = 18))

longitudinal_CFU_w4
ggsave(file.path(plot_dir, "Rebound longitudinal CFU Week 4 HRZE.jpg"),
       longitudinal_CFU_w4, width = 7, height = 5)

# ---- 16S plot (Week 4) ----
rna_w4_16S_plot <- filter(rna_w4_16S, Group %in% c("UNTX", "HRZE 2 WK DOSE", "HRZE 4 WK DOSE"))

longitudinal_16S_w4 <- ggplot(
  rna_w4_16S_plot,
  aes(
    Days.since.Tx.start, Adjusted.16S.Original,
    color = factor(Group2, levels = c("PreRx", "HRZE 2 WK DOSE", "HRZE 4 WK DOSE", "HRZE 8 WK DOSE"))
  )
) +
  geom_jitter(width = 0.0, size = 5, shape = 17) +
  geom_line(data = sum_16S_w4, aes(Days.since.Tx.start, Adjusted.16S.Original, group = Group, col = Group), size = 2) +
  geom_line(data = sum_16S2_w4, aes(Days.since.Tx.start, Adjusted.16S.Original, group = Group2, col = Group2), size = 2) +
  scale_y_log10(limits = c(6e+06, 1e+10), breaks = c(1e+04,1e+05,1e+06,1e+07,1e+8,1e+09,1e+10),
                labels = c(4,5,6,7,8,9,10)) +
  annotate("text", x = 5, y = 0, label = "Treatment",
           hjust = 0, vjust = -27.5, color = "black", size = 5, fontface = "bold") +
  annotate("text", x = 5, y = 0, label = " Recovery Phase",
           hjust = -2, vjust = -27.5, color = "black", size = 5, fontface = "bold") +
  scale_x_continuous(limits = c(0,56), breaks = c(0,14,28,42,56,70,84,98,112),
                     labels = c("0", "2","4","6","8","10","12 weeks","14 weeks","16 weeks")) +
  scale_color_manual(values = c("PreRx" = "brown", "On Treatment" = "brown", "HRZE 4 WK DOSE" = "grey25"),
                     labels = c("PreRx" = "On Treatment", "HRZE 4 WK DOSE" = "Off Treatment")) +
  labs(y = "16S rRNA (log10)", x = "Weeks", colour = "Treatment Status") +
  geom_hline(yintercept = 223304578, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 26, linetype = "dashed", color = "red3") +
  theme_bw() +
  theme(legend.position = "", axis.title = element_text(family = "", color = "black", size = 18)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(axis.text.x = element_text(size = 28, color = "black", face = "bold"),
        axis.text.y = element_text(size = 28, color = "black", face = "bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5, size = 18))

longitudinal_16S_w4
ggsave(file.path(plot_dir, "Rebound longitudinal 16SrRNA Week 4 HRZE.jpg"),
       longitudinal_16S_w4, width = 7, height = 5)

plot_grid(longitudinal_CFU_w4, longitudinal_16S_w4, align = "hv", nrow = 1)

##########################################################################################################################
############################################## STATISTICAL ANALYSIS ######################################################
##########################################################################################################################




##########################################################################################################################
############################################## WILCOXON (GROUP4) ##########################################################
##########################################################################################################################

studysubset_day <- rna_rebound_1 %>%
  mutate(Group4 = paste(Group, Washout.period, sep = "_"))

manual_order <- c(
  "UNTX_0",
  "HRZE 2 WK DOSE_1", "HRZE 2 WK DOSE_5", "HRZE 2 WK DOSE_7", "HRZE 2 WK DOSE_11", "HRZE 2 WK DOSE_14", "HRZE 2 WK DOSE_21", "HRZE 2 WK DOSE_28",
  "HRZE 4 WK DOSE_1", "HRZE 4 WK DOSE_4", "HRZE 4 WK DOSE_7", "HRZE 4 WK DOSE_11", "HRZE 4 WK DOSE_14", "HRZE 4 WK DOSE_21", "HRZE 4 WK DOSE_28",
  "HRZE 8 WK DOSE_1", "HRZE 8 WK DOSE_4", "HRZE 8 WK DOSE_7", "HRZE 8 WK DOSE_11", "HRZE 8 WK DOSE_14", "HRZE 8 WK DOSE_21", "HRZE 8 WK DOSE_28",
  "HRZE.HR 12 WK DOSE_1", "HRZE.HR 12 WK DOSE_4", "HRZE.HR 12 WK DOSE_7", "HRZE.HR 12 WK DOSE_11", "HRZE.HR 12 WK DOSE_13", "HRZE.HR 12 WK DOSE_21", "HRZE.HR 12 WK DOSE_28"
)

studysubset_day$Group4 <- factor(studysubset_day$Group4, levels = manual_order)
studysubset_day <- studysubset_day[order(studysubset_day$Group4), ]

perform_wilcox <- function(var) {
  pw <- pairwise.wilcox.test(studysubset_day[[var]], studysubset_day$Group4, p.adjust.method = "none")
  td <- tidy(pw)
  td$test_type <- "wilcoxon signed rank test"
  td$variable  <- var
  td
}

dataset_3 <- bind_rows(
  perform_wilcox("CFU"),
  perform_wilcox("Adjusted.16S.Original")
)

group1_order <- levels(studysubset_day$Group4)

dataset_3$group1 <- factor(dataset_3$group1, levels = group1_order)
dataset_3$group2 <- factor(dataset_3$group2, levels = group1_order)

wide_dataset <- dataset_3 %>%
  tidyr::spread(key = group1, value = p.value) %>%
  arrange(variable) %>%
  select(group2, all_of(intersect(group1_order, names(.))), variable, test_type)

write.csv(wide_dataset, file.path(stats_dir, "Rebound wilcoxon pvalue.csv"), row.names = FALSE)

##########################################################################################################################
############################################## FIGURE COMBINATIONS #######################################################
##########################################################################################################################

make_img_panel <- function(path, label, hjust_val, vjust_val) {
  im <- readJPEG(path)
  ggplot() +
    annotation_custom(rasterGrob(im), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    annotate("text", x = 0, y = 1, label = paste0(label, ")"),
             hjust = hjust_val, vjust = vjust_val, size = 5, fontface = "bold") +
    theme_void()
}

plot1 <- make_img_panel(file.path(plot_dir, "Rebound longitudinal CFU Week 2 HRZE.jpg"),     "a", 18, -15.5)
plot2 <- make_img_panel(file.path(plot_dir, "Rebound longitudinal 16SrRNA Week 2 HRZE.jpg"), "b", 17, -15.5)
plot3 <- make_img_panel(file.path(plot_dir, "Rebound longitudinal CFU Week 4 HRZE.jpg"),     "c", 18, -15.5)
plot4 <- make_img_panel(file.path(plot_dir, "Rebound longitudinal 16SrRNA Week 4 HRZE.jpg"), "d", 17, -15.5)

rebound_longitudinal <- plot_grid(plot1, plot2, plot3, plot4, align = "hv", nrow = 2)
ggsave(file.path(plot_dir, "Rebound longitudinal.jpg"), rebound_longitudinal, width = 14, height = 10)

##########################################################################################################################
############################################## DURATION COMPARISONS ######################################################
##########################################################################################################################

rna_rebound_Z <- rna_rebound_1 %>%
  filter(Days.since.Tx.start %in% c(0, 11, 12, 25, 26, 39, 53, 54)) %>%
  mutate(
    since.TX.Start = case_when(
      Days.since.Tx.start == 0  ~ "PreRx",
      Days.since.Tx.start == 12 ~ "2 wks HRZE",
      Days.since.Tx.start == 26 ~ "4 wks HRZE",
      Days.since.Tx.start == 25 ~ "2 wks HRZE + 2 wk holiday",
      Days.since.Tx.start == 39 & Treatment.days == 11 ~ "2 wks HRZE + 4 wk holiday",
      Days.since.Tx.start == 39 & Treatment.days == 25 ~ "4 wks HRZE + 2 wk holiday",
      Days.since.Tx.start == 54 ~ "8 wks HRZE",
      Days.since.Tx.start == 53 ~ "4 wks HRZE + 4 wk holiday",
      TRUE ~ NA_character_
    ),
    Tx.Group = case_when(
      Washout.period == 0           ~ "PreRx",
      Washout.period == 1           ~ "Full treatment",
      Washout.period %in% c(28, 14) ~ "Half Tx and Half Holiday ",
      TRUE ~ NA_character_
    )
  )

study_log <- rna_rebound_Z %>%
  mutate(
    logCFU = log10(CFU),
    log16S = log10(Adjusted.16S.Original)
  )

PreRxCFU <- median(study_log$logCFU[study_log$Group == "UNTX"], na.rm = TRUE)
PreRx16S <- median(study_log$log16S[study_log$Group == "UNTX"], na.rm = TRUE)

study_log <- study_log %>%
  mutate(
    PreRxCFU = PreRxCFU,
    PreRx16S = PreRx16S,
    log10CFUreduction = logCFU - PreRxCFU,
    log1016Sreduction = log16S - PreRx16S,
    Gap = log1016Sreduction - log10CFUreduction
  )

summary_by_time_ordered <- study_log %>%
  mutate(
    Timepoint = case_when(
      since.TX.Start == "PreRx"                     ~ "Pre-treatment baseline",
      since.TX.Start == "2 wks HRZE"                ~ "Upon completion of 2-week treatment",
      since.TX.Start == "2 wks HRZE + 2 wk holiday" ~ "Upon completion of 2-week treatment followed by 2-week recovery",
      since.TX.Start == "2 wks HRZE + 4 wk holiday" ~ "Upon completion of 2-week treatment followed by a 4-week recovery",
      since.TX.Start == "4 wks HRZE"                ~ "Upon completion of 4-week treatment",
      since.TX.Start == "4 wks HRZE + 4 wk holiday" ~ "Upon completion of 4-week treatment followed by a 4-week recovery",
      since.TX.Start == "8 wks HRZE"                ~ "Upon completion of 8-week treatment",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Timepoint)) %>%
  group_by(Timepoint) %>%
  summarise(across(where(is.numeric), median, na.rm = TRUE), .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

write.csv(summary_by_time_ordered, file.path(stats_dir, "summary_by_time_ordered.csv"), row.names = FALSE)

long_DF <- study_log %>%
  pivot_longer(cols = c(log10CFUreduction, log1016Sreduction),
               names_to = "PD", values_to = "Reducion") %>%
  mutate(PD = factor(PD,
                     levels = c("log10CFUreduction", "log1016Sreduction"),
                     labels = c("CFU", "16S rRNA")))

# Week 4 comparison plot (days 25,26)
long_DF_4 <- long_DF %>% filter(Days.since.Tx.start %in% c(25, 26))

summary_CFU_day_4 <- long_DF_4 %>%
  group_by(since.TX.Start, PD) %>%
  summarise(Reducion = median(Reducion, na.rm = TRUE), .groups = "drop") %>%
  as.data.frame()

HRZE_CFU_wk4 <- ggplot(long_DF_4, aes(x = since.TX.Start, y = Reducion)) +
  geom_point(aes(col = Tx.Group, shape = factor(PD)), size = 4, position = position_jitter(w = 0.1, h = 0)) +
  geom_point(data = summary_CFU_day_4, size = 10, shape = 95) +
  scale_y_continuous(limits = c(-5.5, 0), breaks = c(-5,-4,-3,-2,-1,0)) +
  scale_color_manual(values = c("Full treatment" = "brown","On Treatment" = "brown", "Half Tx and Half Holiday " = "grey25"),
                     labels = c("PreRx" = "On Treatment", "HRZE 4 WK DOSE" = "Off Treatment")) +
  xlab("") +
  scale_x_discrete(breaks = c("2 wks HRZE + 2 wk holiday", "4 wks HRZE"),
                   labels = c("2 wks HRZE\n+\n2 wk holiday", "4 wks HRZE")) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1) +
  ylab("Reduction from baseline (log10)") + labs(title = "") +
  facet_wrap(~PD, strip.position = "top", nrow = 1) +
  theme_bw() +
  theme(legend.position = "") +
  theme(strip.background = element_rect(fill = "white")) +
  theme(strip.text = element_text(colour = "black", size = 1, face = "bold")) +
  theme(strip.text.x = element_text(size = 15), strip.text.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold")) +
  theme(axis.title.x = element_text(size = 1), axis.title.y = element_text(size = 15)) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.24),
    axis.text.x = element_text(size = 12, color = "black", face = "bold"),
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    axis.line = element_blank(),
    axis.ticks = element_line(size = .24),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.length = unit(2, "pt")
  )

HRZE_CFU_wk4
ggsave(file.path(plot_dir, "CFU vs 16S rRNA vased on Treatment duration WK4.jpg"),
       HRZE_CFU_wk4, width = 6, height = 6)

# Week 8 comparison plot (days 53,54)
long_DF_8 <- long_DF %>% filter(Days.since.Tx.start %in% c(53, 54))

summary_CFU_day_8 <- long_DF_8 %>%
  group_by(since.TX.Start, PD) %>%
  summarise(Reducion = median(Reducion, na.rm = TRUE), .groups = "drop") %>%
  as.data.frame()

HRZE_CFU_wk8 <- ggplot(long_DF_8, aes(x = since.TX.Start, y = Reducion)) +
  geom_point(aes(col = Tx.Group, shape = factor(PD)), size = 4, position = position_jitter(w = 0.1, h = 0)) +
  geom_point(data = summary_CFU_day_8, size = 10, shape = 95) +
  scale_y_continuous(limits = c(-5.5,0), breaks = c(-5,-4,-3,-2,-1,0)) +
  scale_color_manual(values = c("Full treatment" = "brown","On Treatment" = "brown", "Half Tx and Half Holiday " = "grey25"),
                     labels = c("PreRx" = "On Treatment", "HRZE 4 WK DOSE" = "Off Treatment")) +
  xlab("") +
  scale_x_discrete(breaks = c("4 wks HRZE + 4 wk holiday", "8 wks HRZE"),
                   labels = c("4 wks HRZE\n+\n4 wk holiday", "8 wks HRZE")) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1) +
  ylab("Reduction from baseline (log10)") + labs(title = "") +
  facet_wrap(~PD, strip.position = "top", nrow = 1) +
  theme_bw() +
  theme(legend.position = "") +
  theme(strip.background = element_rect(fill = "white")) +
  theme(strip.text = element_text(colour = "black", size = 1, face = "bold")) +
  theme(strip.text.x = element_text(size = 15), strip.text.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold")) +
  theme(axis.title.x = element_text(size = 1), axis.title.y = element_text(size = 15)) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.24),
    axis.text.x = element_text(size = 12, color = "black", face = "bold"),
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    axis.line = element_blank(),
    axis.ticks = element_line(size = .24),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.length = unit(2, "pt")
  )

HRZE_CFU_wk8
ggsave(file.path(plot_dir, "CFU vs 16S rRNA based on Treatment duration WK8.jpg"),
       HRZE_CFU_wk8, width = 6, height = 6)

##########################################################################################################################
############################################## LOG10 SUMMARY TABLE #######################################################
##########################################################################################################################

summary_dataset <- rna_rebound_1 %>%
  mutate(
    log10_CFU = log10(CFU),
    log10_Adjusted_16S_Original = log10(Adjusted.16S.Original)
  ) %>%
  group_by(Group, Washout.period) %>%
  summarise(
    median_log10_CFU = median(log10_CFU, na.rm = TRUE),
    median_log10_Adjusted_16S_Original = median(log10_Adjusted_16S_Original, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_dataset, file.path(stats_dir, "Rebound study log10 summary.csv"), row.names = FALSE)

##########################################################################################################################
############################################## FOLD CHANGE (EOT vs RECOVERY) #############################################
##########################################################################################################################

tp_med <- study_log %>%
  group_by(since.TX.Start) %>%
  summarise(
    CFU_med  = median(CFU, na.rm = TRUE),
    r16S_med = median(Adjusted.16S.Original, na.rm = TRUE),
    .groups  = "drop"
  )


fold_table <- tibble::tibble(
  Comparison = c("2wk HRZE -> +4wk holiday", "4wk HRZE -> +4wk holiday"),
  EOT = c("2 wks HRZE", "4 wks HRZE"),
  REC = c("2 wks HRZE + 4 wk holiday", "4 wks HRZE + 4 wk holiday")
) %>%
  left_join(tp_med, by = c("EOT" = "since.TX.Start")) %>%
  rename(CFU_EOT = CFU_med, r16S_EOT = r16S_med) %>%
  left_join(tp_med, by = c("REC" = "since.TX.Start")) %>%
  rename(CFU_REC = CFU_med, r16S_REC = r16S_med) %>%
  mutate(
    CFU_ratio  = CFU_REC  / CFU_EOT,
    r16S_ratio = r16S_REC / r16S_EOT,
    
    CFU_direction  = ifelse(CFU_ratio  >= 1, "Increase", "Decrease"),
    r16S_direction = ifelse(r16S_ratio >= 1, "Increase", "Decrease"),
    
    CFU_fold  = ifelse(CFU_ratio  >= 1, round(CFU_ratio,  2), round(1 / CFU_ratio,  2)),
    r16S_fold = ifelse(r16S_ratio >= 1, round(r16S_ratio, 2), round(1 / r16S_ratio, 2))
  ) %>%
  select(
    Comparison,
    CFU_EOT, CFU_REC, CFU_direction, CFU_fold,
    r16S_EOT, r16S_REC, r16S_direction, r16S_fold
  )

print(fold_table)


write.csv(
  fold_table,
  file.path(stats_dir, "Rebound_EOT_vs_Recovery_fold_change.csv"),
  row.names = FALSE
)


#############################################
############## The End ######################
#############################################
