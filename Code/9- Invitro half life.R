#######################################
##
## Script name: In vitro H2O2 16S rRNA half-life (interval + overall)
##
## Purpose of script: (1) Build a median 16S rRNA time series for H2O2 kill (day 0 from PreRx (H2O2))
##                    (2) Compute interval half-life between consecutive timepoints
##                    (3) Compute overall half-life from OLS of ln(median 16S) vs day
##                    (4) Save a clean CSV table
##
## Author: Samuel Tabor
##
#######################################

###############################################
####################STARTUP####################
###############################################

options(stringsAsFactors = FALSE)

library(readxl)
library(dplyr)

file_path  <- "../DataRaw/Manuscript Datasets.xlsx"
sheet_name <- "Invitro"

group_day0 <- "PreRx (H2O2)"
group_h2o2 <- "H2O2"
max_day    <- 7

out_dir <- "../DataProcessed/Analysis Data/6- Invitro/"
out_csv <- file.path(out_dir, "Invitro_H2O2_16S_half_life_hours.csv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

###############################################
##################IMPORT DATA##################
###############################################

invitro_data <- read_excel(file_path, sheet = sheet_name)

# Robust column name handling (your file shows Treatment.da)
if (!"Treatment.days" %in% names(invitro_data)) {
  td <- names(invitro_data)[grepl("^Treatment\\.d", names(invitro_data))]
  if (length(td) > 0) names(invitro_data)[names(invitro_data) == td[1]] <- "Treatment.days"
}

###############################################
#################CALCULATIONS##################
###############################################

dat <- invitro_data %>%
  transmute(
    Group                 = as.character(Group),
    Treatment.days        = as.numeric(`Treatment.days`),
    Adjusted.16S.Original = as.numeric(Adjusted.16S.Original)
  ) %>%
  filter(Group %in% c(group_day0, group_h2o2)) %>%
  filter(!is.na(Treatment.days))

series <- bind_rows(
  dat %>%
    filter(Group == group_day0, Treatment.days == 0) %>%
    summarise(
      Treatment.days = 0,
      med_16S = median(Adjusted.16S.Original, na.rm = TRUE),
      .groups = "drop"
    ),
  dat %>%
    filter(Group == group_h2o2, Treatment.days > 0, Treatment.days <= max_day) %>%
    group_by(Treatment.days) %>%
    summarise(med_16S = median(Adjusted.16S.Original, na.rm = TRUE), .groups = "drop")
) %>%
  arrange(Treatment.days) %>%
  mutate(ln_med_16S = ifelse(!is.na(med_16S) & med_16S > 0, log(med_16S), NA_real_)) %>%
  filter(!is.na(ln_med_16S))

interval_tbl <- series %>%
  mutate(day_next = lead(Treatment.days),
         ln_next  = lead(ln_med_16S)) %>%
  filter(!is.na(day_next), !is.na(ln_next)) %>%
  transmute(
    Type = "Interval",
    `Treatment Days` = day_next,
    dt_days = day_next - Treatment.days,
    k_per_day = - (ln_next - ln_med_16S) / dt_days,
    `Half-Life (hours)` = (log(2) / k_per_day) * 24
  ) %>%
  select(Type, `Treatment Days`, `Half-Life (hours)`)

ols_fit <- lm(ln_med_16S ~ Treatment.days, data = series)
slope_per_day <- unname(coef(ols_fit)[["Treatment.days"]])

overall_tbl <- data.frame(
  Type = "Overall",
  `Treatment Days` = NA_real_,                       # keep numeric so bind_rows works
  `Half-Life (hours)` = (log(2) / (-slope_per_day)) * 24,
  check.names = FALSE
)

half_life_table <- bind_rows(interval_tbl, overall_tbl) %>%
  mutate(
    `Half-Life (hours)` = round(`Half-Life (hours)`, 1),
    `Treatment Days` = ifelse(is.na(`Treatment Days`), "", as.character(`Treatment Days`))
  )

write.csv(half_life_table, out_csv, row.names = FALSE)

half_life_table
