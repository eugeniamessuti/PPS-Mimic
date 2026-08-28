# ============================================================
# Create analysis dataset from MIMIC-IV CSV files
# ============================================================

library(dplyr)
library(readr)
library(lubridate)

setwd("C:/Users/messuti/Documents/MIMIC-IV")

# Read files
patients <- read.csv("patients.csv")
admissions <- read.csv("admissions.csv")
micro <- read.csv("microbiologyevents.csv")


# Keep only needed columns
patients <- patients %>%
  select(subject_id, gender, anchor_age, anchor_year)

admissions <- admissions %>%
  select(subject_id, hadm_id, admittime, dischtime) %>%
  mutate(
    admittime = ymd_hms(admittime),
    dischtime = ymd_hms(dischtime)
  )

micro <- micro %>%
  select(subject_id, hadm_id, charttime, org_name) %>%
  mutate(
    charttime = ymd_hms(charttime)
  )


# Combine admissions and patients
# Keep only hospitalizations longer than 48 hours
my_mimic <- admissions %>%
  left_join(patients, by = "subject_id") %>%
  mutate(
    age_at_admission = anchor_age + year(admittime) - anchor_year
  ) %>%
  filter(
    dischtime > admittime + hours(48)
  )


# Find first positive microbiology event
# occurring MORE THAN 48 hours after admission
hai <- micro %>%
  inner_join(
    my_mimic %>% select(subject_id, hadm_id, admittime),
    by = c("subject_id", "hadm_id")
  ) %>%
  filter(
    !is.na(org_name),
    as.numeric(difftime(charttime, admittime, units = "hours")) > 48
  ) %>%
  group_by(subject_id, hadm_id) %>%
  slice_min(
    charttime,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  select(
    subject_id,
    hadm_id,
    hai_time = charttime
  )


# Add HAI information to main data frame
my_mimic <- my_mimic %>%
  left_join(
    hai,
    by = c("subject_id", "hadm_id")
  ) %>%
  mutate(
    hai = ifelse(is.na(hai_time), 0, 1)
  ) %>%
  select(
    subject_id,
    hadm_id,
    admittime,
    dischtime,
    gender,
    age_at_admission,
    hai,
    hai_time
  )


# Add LOS and time to HAI in days
my_mimic$los <- as.numeric(
  difftime(
    my_mimic$dischtime,
    my_mimic$admittime,
    units = "days"
  )
)

my_mimic$time_to_hai <- as.numeric(
  difftime(
    my_mimic$hai_time,
    my_mimic$admittime,
    units = "days"
  )
)


# Censor LOS at day 60
# 1 = observed, 0 = censored
my_mimic$cens_los <- ifelse(
  my_mimic$los <= 60,
  1,
  0
)

my_mimic$cens_los_time <- pmin(
  my_mimic$los,
  60
)


# Censor HAI at day 60
# 1 = HAI observed within 60 days
# 0 = censored
my_mimic$cens_hai <- ifelse(
  !is.na(my_mimic$time_to_hai) &
    my_mimic$time_to_hai <= 60,
  1,
  0
)

my_mimic$cens_hai_time <- ifelse(
  my_mimic$cens_hai == 1,
  my_mimic$time_to_hai,
  pmin(my_mimic$los, 60)
)


# Save final dataframe
write.csv(
  my_mimic,
  "my_mimic.csv",
  row.names = FALSE
)
