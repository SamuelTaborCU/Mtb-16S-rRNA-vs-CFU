#######################################
##
## Script name: Single drug treated mice
##
## Purpose of script: Generate manuscript plots/tables for single-drug study
##
## Author: Samuel Tabor
##
#######################################

## load the packages we will need:
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
studysubset <- read_excel(file_path, sheet = "Single drug")

out_stats_dir <- "../DataProcessed/Analysis Data/2-Single drugs"
out_plot_dir  <- "../DataProcessed/Plots/2-Single drugs"
dir.create(out_stats_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_plot_dir,  recursive = TRUE, showWarnings = FALSE)

#######################################
## 02 - Data processing
#######################################

group_levels <- c("PreRx", "UNTX", "BDQ25", "BDQ5", "RIF10", "RIF30", "STR", "EMB", "PZA", "INH")

studysubset <- studysubset %>%
  mutate(
    Group = factor(Group, levels = group_levels),
    CFU = as.numeric(CFU)
  ) %>%
  arrange(Group, Treatment.days)

studysubset_day <- studysubset %>%
  filter(Treatment.days %in% c(0, 28))

summary_CFU_day <- studysubset_day %>%
  group_by(Group) %>%
  summarise(CFU = mean(CFU, na.rm = TRUE), .groups = "drop")

summary_16s_day <- studysubset_day %>%
  group_by(Group) %>%
  summarise(Adjusted.16S.Original = mean(Adjusted.16S.Original, na.rm = TRUE), .groups = "drop")

#######################################
## 03 - Change from PreRx (correlation + regression)
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
    CFU_decline     = prerx_means$log10_CFU - log10_CFU,
    Adj16S_decline  = prerx_means$log10_16S - log10_16S,
    Gap             = CFU_decline - Adj16S_decline
  )

cr  <- with(studysubset_day1, cor.test(CFU_decline, Adj16S_decline))
fit <- lm(Adj16S_decline ~ CFU_decline, data = studysubset_day1)

slope     <- coef(fit)[2]
r2        <- summary(fit)$r.squared
pval      <- summary(fit)$coefficients[2, 4]
pval_text <- ifelse(pval < 0.001, "<0.001", round(pval, 6))

slope_ci  <- confint(fit)[2, ]
ci_low    <- slope_ci[1]
ci_up     <- slope_ci[2]

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
  "../DataProcessed/Analysis Data/2-Single drugs/single drugs change from PreRx correlation and regression.csv",
  row.names = FALSE
)

#######################################
## 04 - Log reduction summary table
#######################################

B_min <- 0.1   # for CFU = 0

study_log <- studysubset_day %>%
  select(Group, CFU, Adjusted.16S.Original) %>%
  rename(
    B = CFU,
    S = Adjusted.16S.Original
  ) %>%
  mutate(
    B = ifelse(B == 0, B_min, B)
  ) %>%
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

group_order <- c("BDQ5", "BDQ25", "RIF10", "RIF30", "STR", "INH", "PZA", "EMB")

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
    Fold_SN_over_BN = 10^Gap_log10,
    Fold_Reduction  = 10^abs(Gap_log10),
    Fold_Reduction_Text = dplyr::case_when(
      Gap_log10 > 0 ~ paste0(round(Fold_Reduction, 1), "-times greater reduction in CFU"),
      Gap_log10 < 0 ~ paste0(round(Fold_Reduction, 1), "-times smaller reduction in CFU"),
      TRUE          ~ "No difference in reduction"
    ),
    Group = factor(Group, levels = group_order)
  ) %>%
  arrange(Group) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 2))
  )

write.csv(
  summary_wide,
  "../DataProcessed/Analysis Data/2-Single drugs/Single drugs summary log reduction table.csv",
  row.names = FALSE,
  quote = TRUE
)

# long_DF needed by the existing plot code below (no changes to plot code)
long_DF <- study_log %>%
  mutate(Change_CFU = dlogB, Change_Adj16S = dlogS) %>%
  pivot_longer(
    cols = c(Change_CFU, Change_Adj16S),
    names_to = "PD",
    values_to = "Change"
  ) %>%
  mutate(PD = recode(PD, Change_CFU = "CFU", Change_Adj16S = "Adj16S"))

#######################################
## 05 - Plot: CFU vs 16S log reduction from PreRx
#######################################

# Summarize mean reduction per group
summary_CFU_dayx <- long_DF %>%
  filter(Group != "PreRx") %>%
  group_by(Group, PD) %>%
  summarise(
    Change = mean(Change, na.rm = TRUE),
    .groups = "drop"
  )
long_DF1<- long_DF %>%
  filter(Group != "PreRx")
# Plot
HRZE_CFU <- ggplot(long_DF1, aes(x = PD, y = Change)) +
  geom_point(aes(colour = Group, shape = factor(PD)),
             size = 3, position = position_jitter(width = 0.1, height = 0)) +
  geom_point(data = summary_CFU_dayx, aes(x = PD, y = Change),
             size = 10, shape = 95) +
  scale_y_continuous(limits = c(-4, 1), breaks = -4:1) +
  scale_x_discrete(labels = c(CFU = "CFU", Adj16S = "16S")) +
  scale_colour_manual(values = c(
    "PZA" = "chocolate1", "EMB" = "aquamarine4", "STR" = "brown1",
    "RIF10" = "red", "RIF30" = "red4", "INH" = "blue",
    "BDQ25" = "darkorchid4", "BDQ5" = "darkorchid"
  )) +
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed", linewidth = 1) +
  labs(x = NULL, y = "Change from baseline (log10)", title = NULL) +
  facet_wrap(~ Group, strip.position = "top", nrow = 1) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(colour = "black", size = 25, face = "bold"),
    axis.text.x = element_text(size = 25, colour = "black", face = "bold"),
    axis.text.y = element_text(size = 25, colour = "black", face = "bold"),
    axis.title.x = element_text(size = 1),
    axis.title.y = element_text(size = 25, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.24),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.24),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

HRZE_CFU
ggsave(
  filename = file.path(out_plot_dir, "Single Drugs CFU vs 16S rRNA log Reduction from PreRx.eps"),
  plot = HRZE_CFU,
  device = cairo_ps,
  width = 17,
  height = 6,
  units = "in",
  bg = "white",
  fallback_resolution = 600
)
ggsave(
  filename = file.path(out_plot_dir, "Single Drugs CFU vs 16S rRNA log Reduction from PreRx.tiff"),
  plot = HRZE_CFU,
  device = "tiff",
  width = 17,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

#######################################
## 06 - Pairwise Wilcoxon p-values
#######################################

pairwise_wilcox_tidy <- function(x, g) {
  pw <- pairwise.wilcox.test(x, g, p.adjust.method = "none", paired = FALSE, exact = FALSE)
  mat <- pw$p.value
  if (is.null(mat)) return(tibble(group1 = character(), group2 = character(), p.value = numeric()))
  tibble(
    group1 = rep(colnames(mat), each = nrow(mat)),
    group2 = rep(rownames(mat), times = ncol(mat)),
    p.value = as.numeric(mat)
  ) %>%
    filter(!is.na(p.value))
}

res_cfu  <- pairwise_wilcox_tidy(studysubset_day$log10_CFU, studysubset_day$Group)  %>% mutate(variable = "Log10 CFU")
res_16s  <- pairwise_wilcox_tidy(studysubset_day$log10_16S, studysubset_day$Group)  %>% mutate(variable = "Log10 16S rRNA")

dataset_3 <- bind_rows(res_cfu, res_16s) %>%
  mutate(
    test_type = "Wilcoxon rank-sum test",
    Treatment.days = {
      td <- setdiff(sort(unique(studysubset_day$Treatment.days)), 0); if (length(td)) td[1] else NA_character_
    }
  )

group1_order   <- unique(studysubset_day$Group)
variable_order <- c("Log10 CFU", "Log10 16S rRNA")

dataset_3 <- dataset_3 %>%
  mutate(
    group1 = factor(group1, levels = group1_order),
    group2 = factor(group2, levels = group1_order),
    variable = factor(variable, levels = variable_order),
    Treatment.days = as.character(Treatment.days)
  )

wide_dataset <- dataset_3 %>%
  tidyr::pivot_wider(names_from = group1, values_from = p.value) %>%
  arrange(variable) %>%
  select(group2, all_of(intersect(group1_order, names(.))), variable, Treatment.days, test_type) %>%
  mutate(across(where(is.numeric), ~ signif(.x, 1)))

write.csv(
  wide_dataset,
  "../DataProcessed/Analysis Data/2-Single drugs/Single drugs wilcoxon pvalues.csv",
  row.names = FALSE, na = ""
)

#######################################
## 07 - Ranking plot
#######################################

merged_df <- merge(summary_16s_day, summary_CFU_day, by = "Group")
merged_df$log10_CFU <- log10(merged_df$CFU)
merged_df$log10_16S <- log10(merged_df$Adjusted.16S.Original)

preRx_values <- merged_df %>%
  filter(Group == "PreRx") %>%
  summarize(PreRx_log10_CFU = first(log10_CFU),
            PreRx_log10_16S = first(log10_16S))

merged_df <- merged_df %>%
  mutate(Change_log10_CFU =  preRx_values$PreRx_log10_CFU -log10_CFU,
         Change_log10_16S = preRx_values$PreRx_log10_16S - log10_16S) %>%
  filter(Group != "PreRx")

df_long <- merged_df %>%
  pivot_longer(cols = c(Change_log10_CFU, Change_log10_16S), names_to = "Measurement", values_to = "Value") %>%
  group_by(Measurement) %>%
  mutate(Rank = dense_rank(desc(Value)))

write.csv(df_long, file = "../DataProcessed/Analysis Data/2-Single drugs/single drug ranking.csv", row.names = FALSE)

group_colors <- c("PZA"="chocolate1",
                  "EMB"="aquamarine4","STR"="brown",
                  "RIF"="darkgoldenrod","INH"="blue",
                  "BDQ25"="darkorchid4","RIF10"="red","RIF30"="red4",
                  "BDQ5"="darkorchid")

ranking <- ggplot(df_long, aes(x = Measurement, y = Rank, group = Group, color = Group)) +
  geom_bump(size = 4) +
  geom_label(aes(label =  Group, y = Rank, fill = Group, colour = Group),
             hjust = 0.5, size = 10, fontface = "bold",
             label.padding = unit(0.15, "lines"), fill = "white") +
  scale_colour_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_y_reverse( breaks = c(1,2,3,4,5,6,7,8),labels =c(1,2,3,4,5,6,7,8)) +
  theme_classic() +
  theme(legend.position = "", axis.title = element_text(size=18)) +
  theme(axis.text.x = element_text(size=21, color = "black", face = "bold"),
        axis.text.y = element_text(size=21, color = "black", face = "bold"))+
  labs(title = " ",
       x = "",
       y = "Ranking",
       color = "Drug") +
  scale_x_discrete(limits = c("Change_log10_CFU", "Change_log10_16S"), labels = c("CFU", "16S rRNA")) +
  guides(color = guide_legend(override.aes = list(fill = NA, label = " ")))

print(ranking)
ggsave(
  filename = file.path(out_plot_dir, "Single Drugs ranking.eps"),
  plot = ranking,
  device = cairo_ps,
  width = 7,
  height = 6,
  units = "in",
  bg = "white",
  fallback_resolution = 600
)
ggsave(
  filename = file.path(out_plot_dir, "Single Drugs ranking.tiff"),
  plot = ranking,
  device = "tiff",
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

#######################################
## 08 - Scatter: change from PreRx (Panel C)
#######################################

df_full <- studysubset_day %>%
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

PreRxCFU_log <- mean(df_full$CFU [df_full$Group == "PreRx"], na.rm = TRUE)
PreRx16S_log <- mean(df_full$r16S[df_full$Group == "PreRx"], na.rm = TRUE)

df <- df_full %>%
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

fit  <- lm(d16S ~ dCFU, data = df)

slope  <- unname(coef(fit)[2])
ci     <- confint(fit)[2, ]
ci_low <- ci[1]
ci_up  <- ci[2]

pval <- summary(fit)$coefficients[2, 4]
pval_text <- ifelse(pval < 0.001, "<0.001", sprintf("%.6f", pval))

lab <- sprintf(
  "Slope = %.2f; p = %s\n95%% CI (%.2f–%.2f)",
  slope, pval_text, ci_low, ci_up
)

max_abs <- max(abs(c(df$dCFU, df$d16S)), na.rm = TRUE)
x_range <- c(-max_abs, max_abs)
y_range <- c(-max_abs, max_abs)

x_anno  <- x_range[1] + 0.05 * diff(x_range)
y_anno  <- y_range[2] - 0.05 * diff(y_range)

pc <- ggplot(df, aes(x = dCFU, y = d16S, col = Group)) +
  geom_point(size = 2.8, alpha = 0.9) +
  geom_smooth(
    data = df,
    aes(x = dCFU, y = d16S),
    inherit.aes = FALSE,
    method = "lm", se = FALSE,
    linewidth = 1.5, color = "black"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_continuous(limits = x_range) +
  scale_y_continuous(limits = y_range) +
  scale_colour_manual(values = c(
    "PreRx" = "grey50",
    "PZA"   = "chocolate1",
    "EMB"   = "aquamarine4",
    "STR"   = "brown1",
    "RIF"   = "darkgoldenrod",
    "INH"   = "blue",
    "BDQ25" = "darkorchid4",
    "RIF10" = "red",
    "RIF30" = "red4",
    "BDQ5"  = "darkorchid"
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

pc
ggsave(
  filename = file.path(out_plot_dir, "Single_drugs_panel_B_scatter.eps"),
  plot = pc,
  device = cairo_ps,
  width = 7,
  height = 5,
  units = "in",
  bg = "white",
  fallback_resolution = 600
)
ggsave(
  filename = file.path(out_plot_dir, "Single_drugs_panel_B_scatter.tiff"),
  plot = pc,
  device = "tiff",
  width = 7,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

