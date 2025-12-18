#######################################
##
## Script name: Multi-drug regimen 
##
## Purpose of script: Generate manuscript plots/tables for multi-drug regimens
##
## Author: Samuel Tabor
##
#######################################

## load the packages 
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbump)
library(grid)

#######################################
## 01 - Read in data
#######################################

file_path <- "../DataRaw/Manuscript Datasets.xlsx"
studysubset <- read_excel(file_path, sheet = "Multi drugs regimen")

out_stats_dir <- "../DataProcessed/Analysis Data/3- Regimen treated"
out_plot_dir  <- "../DataProcessed/Plots/3- Regimen treated"
dir.create(out_stats_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_plot_dir,  recursive = TRUE, showWarnings = FALSE)

#######################################
## 02 - Data processing
#######################################

group_order <- c("PreRx", "HRZE", "PaMZ", "BPaL", "BPaMZ",
                 "PZM", "BZM", "BZMRb", "BDOS", "BPaOS")

studysubset <- studysubset %>%
  filter(Group %in% group_order) %>%
  mutate(
    Group = factor(Group, levels = group_order),
    CFU = as.numeric(CFU),
    Adjusted.16S.Original = as.numeric(Adjusted.16S.Original)
  ) %>%
  arrange(Group, Treatment.days)

studysubset_day <- studysubset %>% filter(Treatment.days %in% c(0, 28))

#######################################
## 03 - N per Group x Day (wide)
#######################################

sum_studysubset_wide <- studysubset %>%
  group_by(Group, Treatment.days) %>%
  summarise(N = n(), .groups = "drop") %>%
  pivot_wider(id_cols = Group, names_from = Treatment.days, values_from = N) %>%
  as.data.frame()

print(sum_studysubset_wide)

#######################################
## 04 - Per-study summary (Day 0 and 28; log10 first)
#######################################

summary_CFU_day_study <- studysubset_day %>%
  mutate(logCFU = log10(CFU)) %>%
  group_by(Group, Study.name) %>%
  summarise(
    CFU = mean(CFU, na.rm = TRUE),
    logCFU = mean(logCFU, na.rm = TRUE),
    .groups = "drop"
  )

summary_16s_day_study <- studysubset_day %>%
  mutate(logAdjusted.16S.Original = log10(Adjusted.16S.Original)) %>%
  group_by(Group, Study.name) %>%
  summarise(
    Adjusted.16S.Original = mean(Adjusted.16S.Original, na.rm = TRUE),
    logAdjusted.16S.Original = mean(logAdjusted.16S.Original, na.rm = TRUE),
    .groups = "drop"
  )

Summary_day_study <- merge(summary_CFU_day_study, summary_16s_day_study, by = c("Group", "Study.name"))

PreRx_values_study <- Summary_day_study %>%
  filter(Group == "PreRx") %>%
  select(
    Study.name,
    PreRx_CFU = CFU,
    PreRx_logCFU = logCFU,
    PreRx_16s = Adjusted.16S.Original,
    PreRx_log16s = logAdjusted.16S.Original
  )

summary_result_study <- merge(Summary_day_study, PreRx_values_study, by = "Study.name") %>%
  mutate(
    CFU_Decline_from_PreRx        = PreRx_CFU - CFU,
    CFU_log10_Decline_from_PreRx  = round(PreRx_logCFU - logCFU, 2),
    r16S_Decline_from_PreRx       = PreRx_16s - Adjusted.16S.Original,
    r16S_log10_Decline_from_PreRx = round(PreRx_log16s - logAdjusted.16S.Original, 2),
    CFU_fold_decline              = round(PreRx_CFU / CFU, 2),
    r16S_fold_decline             = round(PreRx_16s / Adjusted.16S.Original, 2),
    Group = factor(Group, levels = group_order)
  ) %>%
  arrange(Group)

print(summary_result_study)

write.csv(
  summary_result_study,
  file = file.path(out_stats_dir, "Regimen_per_Study_change_from_PreRx_Summary.csv"),
  row.names = FALSE
)

#######################################
## 05 - Combined summary across studies (Day 0 and 28; log10 first)
#######################################

summary_CFU_day_all <- studysubset_day %>%
  mutate(logCFU = log10(CFU)) %>%
  group_by(Group) %>%
  summarise(
    CFU = mean(CFU, na.rm = TRUE),
    logCFU = mean(logCFU, na.rm = TRUE),
    .groups = "drop"
  )

summary_16s_day_all <- studysubset_day %>%
  mutate(logAdjusted.16S.Original = log10(Adjusted.16S.Original)) %>%
  group_by(Group) %>%
  summarise(
    Adjusted.16S.Original = mean(Adjusted.16S.Original, na.rm = TRUE),
    logAdjusted.16S.Original = mean(logAdjusted.16S.Original, na.rm = TRUE),
    .groups = "drop"
  )

Summary_day_all <- merge(summary_CFU_day_all, summary_16s_day_all, by = "Group")

PreRx_values_all <- Summary_day_all %>%
  filter(Group == "PreRx") %>%
  select(
    PreRx_CFU = CFU,
    PreRx_logCFU = logCFU,
    PreRx_16s = Adjusted.16S.Original,
    PreRx_log16s = logAdjusted.16S.Original
  )

summary_result_all <- Summary_day_all %>%
  mutate(
    CFU_Decline_from_PreRx        = PreRx_values_all$PreRx_CFU - CFU,
    CFU_log10_Decline_from_PreRx  = round(PreRx_values_all$PreRx_logCFU - logCFU, 2),
    r16S_Decline_from_PreRx       = PreRx_values_all$PreRx_16s - Adjusted.16S.Original,
    r16S_log10_Decline_from_PreRx = round(PreRx_values_all$PreRx_log16s - logAdjusted.16S.Original, 2),
    CFU_fold_decline              = round(PreRx_values_all$PreRx_CFU / CFU, 2),
    r16S_fold_decline             = round(PreRx_values_all$PreRx_16s / Adjusted.16S.Original, 2),
    Group = factor(Group, levels = group_order)
  ) %>%
  arrange(Group)

print(summary_result_all)

write.csv(
  summary_result_all,
  file = file.path(out_stats_dir, "Regimen_combined_studies_change_from_PreRx_Summary.csv"),
  row.names = FALSE
)

#######################################
## 06 - Per-study summary (all timepoints; log10 first)
#######################################

summary_CFU_alltime <- studysubset %>%
  mutate(logCFU = log10(CFU)) %>%
  group_by(Group, Study.name, Treatment.days) %>%
  summarise(
    CFU = mean(CFU, na.rm = TRUE),
    logCFU = mean(logCFU, na.rm = TRUE),
    .groups = "drop"
  )

summary_16s_alltime <- studysubset %>%
  mutate(logAdjusted.16S.Original = log10(Adjusted.16S.Original)) %>%
  group_by(Group, Study.name, Treatment.days) %>%
  summarise(
    Adjusted.16S.Original = mean(Adjusted.16S.Original, na.rm = TRUE),
    logAdjusted.16S.Original = mean(logAdjusted.16S.Original, na.rm = TRUE),
    .groups = "drop"
  )

Summary_alltime <- merge(summary_CFU_alltime, summary_16s_alltime, by = c("Group", "Study.name", "Treatment.days"))

PreRx_values_alltime <- Summary_alltime %>%
  filter(Group == "PreRx") %>%
  select(
    Study.name,
    PreRx_CFU = CFU,
    PreRx_logCFU = logCFU,
    PreRx_16s = Adjusted.16S.Original,
    PreRx_log16s = logAdjusted.16S.Original
  )

summary_result_alltime <- merge(Summary_alltime, PreRx_values_alltime, by = "Study.name") %>%
  mutate(
    CFU_Decline_from_PreRx        = PreRx_CFU - CFU,
    CFU_log10_Decline_from_PreRx  = round(PreRx_logCFU - logCFU, 2),
    r16S_Decline_from_PreRx       = PreRx_16s - Adjusted.16S.Original,
    r16S_log10_Decline_from_PreRx = round(PreRx_log16s - logAdjusted.16S.Original, 2),
    CFU_fold_decline              = round(PreRx_CFU / CFU, 2),
    r16S_fold_decline             = round(PreRx_16s / Adjusted.16S.Original, 2),
    Group = factor(Group, levels = group_order)
  ) %>%
  arrange(Group)

print(summary_result_alltime)

write.csv(
  summary_result_alltime,
  file = file.path(out_stats_dir, "Regimen_per_Study_change_from_PreRx__all_timepoints_Summary.csv"),
  row.names = FALSE
)

#######################################
## 07 - Change from PreRx (correlation + regression; Day 0 and 28)
#######################################

studysubset_day <- studysubset_day %>%
  mutate(
    log10_CFU = log10(CFU),
    log10_16S = log10(Adjusted.16S.Original)
  )

prerx_means <- studysubset_day %>%
  filter(Group == "PreRx") %>%
  summarise(across(c(log10_CFU, log10_16S), ~ mean(.x, na.rm = TRUE)))

studysubset_day1 <- studysubset_day %>%
  filter(Group != "PreRx") %>%
  mutate(
    CFU_decline    = prerx_means$log10_CFU - log10_CFU,
    Adj16S_decline = prerx_means$log10_16S - log10_16S
  )

cr <- with(studysubset_day1, cor.test(CFU_decline, Adj16S_decline))
fit <- lm(Adj16S_decline ~ CFU_decline, data = studysubset_day1)

slope     <- coef(fit)[2]
r2        <- summary(fit)$r.squared
pval      <- summary(fit)$coefficients[2, 4]
pval_text <- ifelse(pval < 0.001, "<0.001", round(pval, 6))

slope_ci <- confint(fit)[2, ]
ci_low   <- slope_ci[1]
ci_up    <- slope_ci[2]

combined_results <- data.frame(
  Correlation_Coefficient = unname(cr$estimate),
  Corr_CI_Lower = unname(cr$conf.int[1]),
  Corr_CI_Upper = unname(cr$conf.int[2]),
  Corr_P_value  = unname(cr$p.value),
  Slope = round(slope, 4),
  Slope_CI_Lower = round(unname(ci_low), 4),
  Slope_CI_Upper = round(unname(ci_up), 4),
  Slope_P_value = pval_text,
  R_squared = round(r2, 4),
  row.names = NULL
)

print(combined_results)

write.csv(
  combined_results,
  file = file.path(out_stats_dir, "Regimen change from PreRx correlation and regression.csv"),
  row.names = FALSE
)

#######################################
## 08 - Log reduction + fold change summary (Day 0 and 28)
#######################################

B_min <- 0.1

study_log <- studysubset_day %>%
  select(Group, CFU, Adjusted.16S.Original) %>%
  rename(B = CFU, S = Adjusted.16S.Original) %>%
  mutate(B = ifelse(B == 0, B_min, B)) %>%
  filter(B > 0, S > 0) %>%
  mutate(
    logB = log10(B),
    logS = log10(S)
  )

baseline <- study_log %>%
  filter(Group == "PreRx") %>%
  summarise(
    logB0 = mean(logB, na.rm = TRUE),
    logS0 = mean(logS, na.rm = TRUE)
  )

study_log <- study_log %>%
  mutate(
    dlogB     = logB - baseline$logB0,
    dlogS     = logS - baseline$logS0,
    Gap_log10 = dlogS - dlogB
  )

summary_wide <- study_log %>%
  filter(Group != "PreRx") %>%
  group_by(Group) %>%
  summarise(
    B_change_log10 = mean(dlogB,     na.rm = TRUE),
    S_change_log10 = mean(dlogS,     na.rm = TRUE),
    Gap_log10      = mean(Gap_log10, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Fold_SN_over_BN = round(10^Gap_log10),
    Fold_Reduction  = round(10^abs(Gap_log10)),
    B_change_log10  = round(B_change_log10, 2),
    S_change_log10  = round(S_change_log10, 2),
    Gap_log10       = round(Gap_log10, 2),
    Fold_Reduction_Text = dplyr::case_when(
      Gap_log10 > 0 ~ paste0(Fold_Reduction, "-times greater reduction in CFU"),
      Gap_log10 < 0 ~ paste0(Fold_Reduction, "-times smaller reduction in CFU"),
      TRUE          ~ "No difference in reduction"
    ),
    Group = factor(Group, levels = group_order)
  ) %>%
  arrange(Group)

print(summary_wide)

write.csv(
  summary_wide,
  file = file.path(out_stats_dir, "Regimen log reduction and fold change summary.csv"),
  row.names = FALSE
)

#######################################
## 09 - Per-study reductions (Day 0 and 28; used for gap table)
#######################################

study_log_ps <- studysubset_day[, c("Group", "CFU", "Adjusted.16S.Original", "Study.name")]
study_log_ps$logCFU <- log10(study_log_ps$CFU)
study_log_ps$log16S <- log10(study_log_ps$Adjusted.16S.Original)

PreRxCFU <- study_log_ps %>%
  filter(Group == "PreRx") %>%
  group_by(Study.name) %>%
  summarise(PreRxCFU = mean(logCFU, na.rm = TRUE), .groups = "drop")

PreRx16S <- study_log_ps %>%
  filter(Group == "PreRx") %>%
  group_by(Study.name) %>%
  summarise(PreRx16S = mean(log16S, na.rm = TRUE), .groups = "drop")

study_log_ps <- left_join(study_log_ps, PreRxCFU, by = "Study.name")
study_log_ps <- left_join(study_log_ps, PreRx16S, by = "Study.name")

study_log_ps <- subset(study_log_ps, Group != "PreRx")
study_log_ps$log10CFUreduction <- study_log_ps$logCFU - study_log_ps$PreRxCFU
study_log_ps$log1016Sreduction <- study_log_ps$log16S - study_log_ps$PreRx16S

study_log_ps <- study_log_ps[, c("Group", "Study.name", "PreRxCFU", "logCFU",
                                 "log10CFUreduction", "PreRx16S", "log16S", "log1016Sreduction")]

long_DF <- pivot_longer(study_log_ps,
                        cols = c(log10CFUreduction, log1016Sreduction),
                        names_to = "PD", values_to = "Reducion")

long_DF$PD <- factor(long_DF$PD,
                     levels = c("log10CFUreduction", "log1016Sreduction"),
                     labels = c("CFU", "a16S"))

summary_CFU_dayx <- long_DF %>%
  group_by(Group, PD, Study.name) %>%
  summarise(Reducion = mean(Reducion, na.rm = TRUE), .groups = "drop") %>%
  as.data.frame()

summary_CFU_dayx_wide <- summary_CFU_dayx %>%
  pivot_wider(names_from = PD, values_from = Reducion) %>%
  as.data.frame()

summary_CFU_dayx_wide$Decline_Gap <- summary_CFU_dayx_wide$`CFU` - summary_CFU_dayx_wide$`a16S`

summary_CFU_dayx_wide <- summary_CFU_dayx_wide %>%
  mutate(
    CFU = round(CFU, 2),
    a16S = round(a16S, 2),
    Fold_Change = round(10^abs(Decline_Gap), 1),
    Fold_Change_Description = ifelse(
      Decline_Gap < 0,
      paste0(Fold_Change, "-times greater reduction in CFU"),
      paste0(Fold_Change, "-times smaller reduction in CFU")
    )
  )

print(summary_CFU_dayx_wide)

#######################################
## 10 - Pairwise Wilcoxon p-values (Day 0 and 28)
#######################################

pairwise_wilcox_tidy <- function(x, g) {
  pw <- pairwise.wilcox.test(x, g, p.adjust.method = "none", paired = FALSE, exact = FALSE)
  mat <- pw$p.value
  if (is.null(mat)) return(tibble(group1 = character(), group2 = character(), p.value = numeric()))
  tibble(
    group1 = rep(colnames(mat), each = nrow(mat)),
    group2 = rep(rownames(mat), times = ncol(mat)),
    p.value = as.numeric(mat)
  ) %>% filter(!is.na(p.value))
}

res_cfu <- pairwise_wilcox_tidy(studysubset_day$log10_CFU, studysubset_day$Group) %>% mutate(variable = "Log10 CFU")
res_16s <- pairwise_wilcox_tidy(studysubset_day$log10_16S, studysubset_day$Group) %>% mutate(variable = "Log10 16S rRNA")

dataset_3 <- bind_rows(res_cfu, res_16s) %>%
  mutate(
    test_type = "Wilcoxon rank-sum test",
    Treatment.days = {
      td <- setdiff(sort(unique(studysubset_day$Treatment.days)), 0)
      if (length(td)) as.character(td[1]) else NA_character_
    }
  )

group1_order   <- unique(studysubset_day$Group)
variable_order <- c("Log10 CFU", "Log10 16S rRNA")

dataset_3 <- dataset_3 %>%
  mutate(
    group1 = factor(group1, levels = group1_order),
    group2 = factor(group2, levels = group1_order),
    variable = factor(variable, levels = variable_order)
  )

wide_dataset <- dataset_3 %>%
  pivot_wider(names_from = group1, values_from = p.value) %>%
  arrange(variable) %>%
  select(group2, all_of(intersect(group1_order, names(.))), variable, Treatment.days, test_type)

p_cols <- intersect(group1_order, names(wide_dataset))
wide_dataset <- wide_dataset %>%
  mutate(across(all_of(p_cols), ~ signif(.x, 1)))

op <- options(scipen = 999); on.exit(options(op), add = TRUE)

write.csv(
  wide_dataset,
  file = file.path(out_stats_dir, "Regimen wilcoxon t-test pvalue.csv"),
  row.names = FALSE, na = ""
)

#######################################
## 11 - Panel d: change-from-PreRx scatter (UNCHANGED PLOT CODE)
#######################################

df <- studysubset_day %>%
  rename(
    day      = Treatment.days,
    r16S_raw = Adjusted.16S.Original
  ) %>%
  mutate(
    CFU_raw  = as.numeric(CFU),
    r16S_raw = as.numeric(r16S_raw),
    day      = as.numeric(day)
  ) %>%
  filter(CFU_raw > 0, r16S_raw > 0) %>%
  mutate(
    CFU  = log10(CFU_raw),
    r16S = log10(r16S_raw)
  )

PreRxCFU_log <- mean(df$CFU [df$Group == "PreRx"], na.rm = TRUE)
PreRx16S_log <- mean(df$r16S[df$Group == "PreRx"], na.rm = TRUE)

df_change <- df %>%
  mutate(
    dCFU = CFU  - PreRxCFU_log,
    d16S = r16S - PreRx16S_log
  ) %>%
  filter(Group != "PreRx")

my_theme <- theme_bw() +
  theme(
    legend.position = "top",
    axis.title = element_text(family = "", color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(
      size  = 18, color = "black", face = "bold",
      angle = 0, hjust = 0.5
    ),
    axis.text.y = element_text(size = 18, color = "black", face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )

fit  <- lm(d16S ~ dCFU, data = df_change)

slope  <- unname(coef(fit)[2])
ci     <- confint(fit)[2, ]
ci_low <- ci[1]
ci_up  <- ci[2]

pval <- summary(fit)$coefficients[2, 4]
pval_text <- ifelse(pval < 0.001, "<0.001", sprintf("%.2f", pval))

lab <- sprintf(
  "Slope = %.2f; p = %s\n95%% CI (%.2f–%.2f)",
  slope, pval_text, ci_low, ci_up
)

max_abs <- 8
x_range <- c(-max_abs, max_abs)
y_range <- c(-max_abs, max_abs)
break_axis <- c(-8, -4, 0, 4, 8)

x_anno  <- x_range[1] + 0.05 * diff(x_range)
y_anno  <- y_range[2] - 0.05 * diff(y_range)

pd <- ggplot(df_change, aes(x = dCFU, y = d16S, col = Group)) +
  geom_point(size = 2.8, alpha = 0.9) +
  geom_smooth(
    data = df_change,
    aes(x = dCFU, y = d16S),
    inherit.aes = FALSE,
    method = "lm", se = FALSE, linewidth = 1.5, color = "black"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_continuous(limits = x_range, breaks = break_axis) +
  scale_y_continuous(limits = y_range, breaks = break_axis) +
  scale_colour_manual(values = c(
    "PreRx"  = "grey50",
    "HRZE"   = "#E4B811",
    "PaMZ"   = "#00A9FF",
    "BDOS"   = "forestgreen",
    "BPaOS"  = "#C77CFF",
    "BPaL"   = "aquamarine4",
    "PZM"    = "black",
    "BZM"    = "gold4",
    "BZMRb"  = "blue",
    "BPaMZ"  = "red",
    "BDO286" = "pink",
    "BDO830" = "blue4",
    "BD286S" = "green",
    "BPa286S"= "orange",
    "BD830S" = "brown1",
    "BPa830S" = "deepskyblue",
    "BPaO286" = "coral3",
    "BOS286" = "bisque",
    "DOS286" = "azure2",
    "PaOS286"= "firebrick"
  )) +
  annotate(
    "text", x = x_anno, y = y_anno, label = lab,
    hjust = 0, vjust = 1, fontface = "bold", size = 7.5
  ) +
  labs(
    x = "Change in CFU from baseline (log10)",
    y = "Change in 16S rRNA from baseline (log10)"
  ) +
  my_theme +
  theme(
    legend.position = "right",
    legend.title    = element_blank(),
    legend.text     = element_text(size = 18, face = "bold"),
    legend.key.size = unit(2, "lines")
  )

pd

ggsave(
  file.path(out_plot_dir, "Regimen_CFU vs_16S_28_panel_d_scatter_without_PreRx.jpg"),
  pd, width = 7, height = 5, units = "in", dpi = 300
)

#######################################
## 12 - Ranking 
#######################################

summary_CFU_day_rank <- studysubset_day %>%
  group_by(Group) %>% summarise(CFU = mean(CFU, na.rm = TRUE), .groups = "drop") %>%
  as.data.frame()

summary_16s_day_rank <- studysubset_day %>%
  group_by(Group) %>% summarise(Adjusted.16S.Original = mean(Adjusted.16S.Original, na.rm = TRUE), .groups = "drop") %>%
  as.data.frame()

merged_df <- merge(summary_16s_day_rank, summary_CFU_day_rank, by = "Group")
merged_df$log10_CFU <- log10(merged_df$CFU)
merged_df$log10_16S <- log10(merged_df$Adjusted.16S.Original)

preRx_values <- merged_df %>%
  filter(Group == "PreRx") %>%
  summarize(
    PreRx_log10_CFU = first(log10_CFU),
    PreRx_log10_16S = first(log10_16S)
  )

merged_df <- merged_df %>%
  mutate(
    Change_log10_CFU = preRx_values$PreRx_log10_CFU - log10_CFU,
    Change_log10_16S = preRx_values$PreRx_log10_16S - log10_16S
  )

write.csv(
  merged_df,
  file = file.path(out_stats_dir, "regimen PD markers Summary.csv"),
  row.names = FALSE
)

merged_df <- merged_df %>% filter(Group != "PreRx")

df_long <- merged_df %>%
  pivot_longer(cols = c(Change_log10_CFU, Change_log10_16S),
               names_to = "Measurement", values_to = "Value") %>%
  group_by(Measurement) %>%
  mutate(Rank = dense_rank(desc(Value)))

write.csv(df_long, file = file.path(out_stats_dir, "Regimen ranking.csv"), row.names = FALSE)

group_colors <- c(
  "PreRx"="grey50","HRZE"= "#E4B811","PaMZ"= "#00A9FF","BDOS"="forestgreen",
  "BPaOS"= "#C77CFF","BPaL"="aquamarine4","PZM"="black",
  "BZM"="gold4","BZMRb"="blue","BPaMZ"="red","BDO286" = "pink",
  "BDO830" = "blue4","BD286S" = "green","BPa286S" = "orange",
  "BD830S" = "brown1","BPa830S" = "deepskyblue","BPaO286" = "coral3",
  "BOS286" = "bisque","DOS286" = "azure2","PaOS286" = "firebrick"
)

ranking <- ggplot(df_long, aes(x = Measurement, y = Rank, group = Group, color = Group)) +
  geom_bump(size = 4) +
  geom_label(aes(label =  Group, y = Rank, fill = Group, colour = Group),
             hjust = 0.5, size = 10, fontface = "bold",
             label.padding = unit(0.15, "lines"), fill = "white") +
  scale_colour_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_y_reverse(breaks = c(1,2,3,4,5,6,7,8,9), labels = c(1,2,3,4,5,6,7,8,9)) +
  theme_classic() +
  theme(legend.position = "", axis.title = element_text(size=18)) +
  theme(axis.text.x = element_text(size=21, color = "black", face = "bold"),
        axis.text.y = element_text(size=21, color = "black", face = "bold"))+
  labs(title = " ",
       x = "",
       y = "Ranking",
       color = "Drug") +
  scale_x_discrete(limits = c("Change_log10_CFU", "Change_log10_16S"),
                   labels = c("CFU", "16S rRNA")) +
  guides(color = guide_legend(override.aes = list(fill = NA, label = " ")))

print(ranking)
ggsave(file.path(out_plot_dir, "Regimen ranking.jpg"), width = 7, height = 6)

#######################################
## 13 - Correlation table (all timepoints; without PreRx)
#######################################

format_p_first_nonzero <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p == 0) return("0")
  dp <- ceiling(-log10(p))
  dp <- max(dp, 1)
  fmt <- paste0("%.", dp, "f")
  sprintf(fmt, round(p, dp))
}

studysubset_dayx <- studysubset %>%
  filter(Group != "PreRx") %>%
  mutate(
    logCFU = log10(CFU),
    log16S = log10(Adjusted.16S.Original)
  ) %>%
  filter(is.finite(logCFU), is.finite(log16S))

correlation_results_all <- studysubset_dayx %>%
  group_by(Treatment.days) %>%
  group_modify(~{
    df <- .x
    n_obs <- nrow(df)
    nC <- dplyr::n_distinct(df$logCFU)
    nS <- dplyr::n_distinct(df$log16S)
    
    if (n_obs >= 3 && nC > 1 && nS > 1) {
      
      res_cor <- cor.test(df$logCFU, df$log16S, method = "pearson")
      corr_num <- as.numeric(unname(res_cor$estimate))
      ci_low   <- as.numeric(res_cor$conf.int[1])
      ci_high  <- as.numeric(res_cor$conf.int[2])
      p_raw    <- as.numeric(res_cor$p.value)
      
      fit <- lm(log16S ~ logCFU, data = df)
      slope <- coef(fit)[2]
      slope_ci <- confint(fit)[2, ]
      slope_low <- slope_ci[1]
      slope_up  <- slope_ci[2]
      r2 <- summary(fit)$r.squared
      
      tibble(
        n = n_obs,
        Corr    = round(corr_num, 2),
        CI_L    = round(ci_low, 2),
        CI_U    = round(ci_high, 2),
        p       = format_p_first_nonzero(p_raw),
        Slope   = round(slope, 2),
        Slope_L = round(slope_low, 2),
        Slope_U = round(slope_up, 2),
        R2      = round(r2, 2)
      )
      
    } else {
      
      tibble(
        n = n_obs,
        Corr = NA_real_,
        CI_L = NA_real_,
        CI_U = NA_real_,
        p    = NA_character_,
        Slope   = NA_real_,
        Slope_L = NA_real_,
        Slope_U = NA_real_,
        R2      = NA_real_
      )
    }
  }) %>%
  ungroup()

print(correlation_results_all)

write.csv(
  correlation_results_all,
  file = file.path(out_stats_dir, "Regimen_CFU vs_16S_Corr without PreRX all timepoints.csv"),
  row.names = FALSE
)

