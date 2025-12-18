#######################################
##
## Script name: Untreated mice
##
## Purpose of script: Generate manuscript plots for untreated mice:
##                    (A) time course of CFU, 16S rRNA, and 16S/CFU (log10)
##                    (B) 16S rRNA vs CFU (log10) with regression
##
## Author: Samuel Tabor
##
#######################################

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(grid)

# Input and Output
file_path <- "../DataRaw/Manuscript Datasets.xlsx"
out_stats_dir <- "../DataProcessed/Analysis Data/1-Untreated mice"
out_plot_dir  <- "../DataProcessed/Plots/1-Untreated mice"
dir.create(out_stats_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_plot_dir,  recursive = TRUE, showWarnings = FALSE)

# Data
STD_2 <- read_excel(file_path, sheet = "Untreated mice")

df <- STD_2 %>%
  rename(
    day = Day.of.Infection.at.tissue.collection,
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

# Summary for Panel A
by_day <- df %>%
  group_by(day) %>%
  summarise(
    CFU_mean  = mean(CFU,  na.rm = TRUE),
    r16S_mean = mean(r16S, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(RNA_CFU = r16S_mean - CFU_mean)

long_day <- by_day %>%
  pivot_longer(
    cols = c(CFU_mean, r16S_mean, RNA_CFU),
    names_to = "measure",
    values_to = "value"
  ) %>%
  mutate(
    measure = recode(measure,
                     CFU_mean  = "CFU",
                     r16S_mean = "16S rRNA",
                     RNA_CFU   = "16SRNA/CFU"),
    measure = factor(measure, levels = c("CFU", "16S rRNA", "16SRNA/CFU"))
  )

# ---- Theme for Panels ----
my_theme <- theme_bw() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 16, face = "bold"),
    legend.key.size = unit(1.2, "cm"),
    axis.title = element_text(family = "", color = "black", size = 20, face = "bold"),
    axis.text.x = element_text(size = 18, color = "black", face = "bold", angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 18, color = "black", face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )

###########################################################
## Panel A
###########################################################
pA <- ggplot(long_day,
             aes(x = day, y = value,
                 color = measure, linetype = measure, shape = measure)) +
  geom_point(size = 5,
             position = position_jitter(w = 0.1, h = 0.01)) +
  geom_line(linewidth = 1.5) +
  
  scale_y_continuous(
    limits = c(1, 9),
    breaks = 1:9,
    name = "CFU or 16S rRNA Count (log10)",
    sec.axis = sec_axis(~ .,  breaks = 1:9,
                        name = "16S rRNA / CFU (log10)")
  ) +
  
  scale_x_continuous(
    limits = c(0, 60),
    breaks = c(0, 7, 11, 19, 28, 56)
  ) +
  
  scale_color_manual(values = c("CFU" = "black",
                                "16S rRNA" = "red3",
                                "16SRNA/CFU" = "dodgerblue4")) +
  scale_linetype_manual(values = c("CFU" = "solid",
                                   "16S rRNA" = "dashed",
                                   "16SRNA/CFU" = "dashed")) +
  scale_shape_manual(values = c("CFU" = 16,
                                "16S rRNA" = 17,
                                "16SRNA/CFU" = 23)) +
  
  labs(x = "Days",
       color = NULL, linetype = NULL, shape = NULL) +
  my_theme

###########################################################
## Panel B 
###########################################################
fit  <- lm(r16S ~ CFU, data = df)

slope <- unname(coef(fit)[2])
ci     <- confint(fit)[2, ]
ci_low <- ci[1]
ci_up  <- ci[2]

pval  <- summary(fit)$coefficients[2, 4]
pval_text <- ifelse(pval < 0.001, "<0.001", sprintf("%.6f", pval))

lab <- sprintf("Slope = %.2f; p = %s\n95%% CI (%.2f–%.2f)",
               slope, pval_text, ci_low, ci_up)

pB <- ggplot(df, aes(x = CFU, y = r16S)) +
  geom_point(size = 5.5, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_continuous(limits = c(3, 9), breaks = 3:9) +
  scale_y_continuous(limits = c(3, 9), breaks = 3:9) +
  annotate("text", x = 4.2, y = 5.0, label = lab,
           hjust = 0, vjust = 1, fontface = "bold", size = 8) +
  labs(x = "CFU (log10)", y = "16S rRNA (log10)") +
  my_theme

# Stats output (matches the plotted data)
correlation_UTD <- cor.test(df$CFU, df$r16S)
r2 <- summary(fit)$r.squared

final_output <- data.frame(
  Correlation_Coefficient = unname(correlation_UTD$estimate),
  CI_Lower = unname(correlation_UTD$conf.int[1]),
  CI_Upper = unname(correlation_UTD$conf.int[2]),
  Correlation_P_value = unname(correlation_UTD$p.value),
  Model = "log10_16S ~ log10_CFU",
  Slope = round(slope, 4),
  Slope_CI_Lower = round(unname(ci_low), 4),
  Slope_CI_Upper = round(unname(ci_up), 4),
  R_squared = round(r2, 4),
  Slope_P_value = pval_text,
  row.names = NULL
)


write.csv(final_output,
          file.path(out_stats_dir, "Combined_Correlation_Regression_Results.csv"),
          row.names = FALSE)

write.csv(
  data.frame(Model = "r16S ~ CFU", Slope = round(slope, 4), P_value = pval_text),
  file.path(out_stats_dir, "16S_vs_CFU_slope_pval.csv"),
  row.names = FALSE
)

ggsave(file.path(out_plot_dir, "Untreated_mice_panel_A_timecourse.jpg"),
       pA, width = 6, height = 5, units = "in", dpi = 300)

ggsave(file.path(out_plot_dir, "Untreated_mice_panel_B_scatter.jpg"),
       pB, width = 6, height = 5, units = "in", dpi = 300)

pA
pB
print(final_output)
