###############################################################################
# calculate the standardised absolute risk (AR) in the cancer and noncancer cohorts

# three options:

# 1. age- and sex-standardised within cohort x vaccine arm
# 2. vaccine arm-, age- and sex-standardised within cohort
# 3. age- and sex-standardised within cohort

###############################################################################
library(tidyverse)
library(glue)

###############################################################################
# setup

# read outcomes
outcomes <- readr::read_rds(
  here::here("analysis", "lib", "outcomes.rds"))

# redaction functions
source(here::here("analysis", "functions", "redaction_functions.R"))

# create output directory
fs::dir_create(here::here("output", "tte", "tables"))

###############################################################################
# read data
# time to event data split into 6 comparison periods
data_tte <- bind_rows(
  lapply(
    outcomes,
    function(x)
      readr::read_rds(
        here::here("output", "tte", "data", glue("data_tte_{x}.rds"))
      ) %>%
      mutate(outcome = x)
  )
) %>%
  mutate(
    cohort = if_else(
      cancer_subgroup %in% "noncancer", 
      cancer_subgroup, 
      "cancer"
      )
    ) %>%
  mutate(across(arm, ~if_else(.x %in% "unvax", as.character(.x), "vax"))) %>%
  select(patient_id, cohort, arm, k, outcome, status)

# define age bands
age_bands <- c(-Inf, seq(25,90,5), Inf)
lower_bands <- as.character(age_bands[-length(age_bands)])
lower_bands[1] <- "18"
upper_bands <- as.character(age_bands[-1] - 1)
upper_bands[length(upper_bands)] <- "120"
age_band_labs <- str_c(lower_bands, upper_bands, sep = "-")

# read age and sex data
data_all <- readr::read_rds(
  here::here("output", "data", "data_all.rds")) %>%
  select(patient_id, age, sex) %>%
  mutate(across(age, ~cut(.x, breaks = age_bands, labels = age_band_labs)))

###############################################################################
### 1. age- and sex-standardised AR within cohort x vaccine arm
# AR grouped by arm, age, sex
data_counts_arm_group <- data_tte %>%
  left_join(data_all, by = "patient_id") %>% 
  group_by(cohort, arm, outcome, k, age, sex) %>%
  summarise(
    n = n(),
    events = sum(status),
    .groups = "keep"
  ) %>%
  ungroup() %>%
  mutate(ar = events/n) 

# within each outcome, arm and comparison period, 
# get proportion of general population in each sex:age_band
weights_arm_group <- data_counts_arm_group %>%
  filter(cohort == "noncancer") %>%
  group_by(outcome, arm, k) %>%
  mutate(prop = n/sum(n)) %>%
  # mutate(total = sum(n)) %>%
  ungroup() %>%
  select(arm, outcome, k, age, sex, weight = prop)

# calculate AR in both arms, 
# weighted by sex:age_band distribution in noncancer pop
data_ar_arm_group <- data_counts_arm_group %>%
  left_join(
    weights_arm_group, by = c("arm", "outcome", "k", "age", "sex")
  ) %>%
  group_by(cohort, arm, outcome, k) %>%
  summarise(
    events = sum(events),
    # calculate crude AR 
    # (this is just to make sure that it is the same as the weighted average in 
    # the noncancer cohort)
    ar_crude = sum(events)/sum(n),
    # calculate the age- and sex-standardised AR
    ar_weighted = weighted.mean(x = ar, w = weight),
    .groups = "keep"
  ) %>%
  ungroup() %>%
  # redact if <=5 events
  # if redact one arm, redact other too to avoid potential for back calculation
  # from previously released results
  mutate(redact = events <= 5) %>%
  group_by(cohort, outcome, k) %>%
  mutate(redact = as.logical(max(redact))) %>%
  ungroup() %>%
  mutate(across(c("events", "ar_crude", "ar_weighted"), 
                ~if_else(redact, NA_real_, as.numeric(.x)))) %>%
  select(-events, -redact)

readr::write_csv(
  data_ar_arm_group,
  here::here("output", "tte", "tables", "ar_arm_group.csv")
  )

### 2. arm-, age- and sex-standardised AR within cohort
# within each outcome and comparison period, 
# get proportion of general population in each arm:sex:age_band
weights_arm_weight <- data_counts_arm_group %>%
  filter(cohort == "noncancer") %>%
  group_by(outcome, k) %>%
  mutate(prop = n/sum(n)) %>%
  # mutate(total = sum(n)) %>%
  ungroup() %>%
  select(arm, outcome, k, age, sex, weight = prop)

# calculate AR, 
# weighted by arm:sex:age_band distribution in noncancer pop
data_ar_arm_weight <- data_counts_arm_group %>%
  left_join(
    weights_arm_weight, by = c("arm", "outcome", "k", "age", "sex")
  ) %>%
  group_by(cohort, outcome, k) %>%
  summarise(
    events = sum(events),
    ar_crude = sum(events)/sum(n),
    ar_weighted = weighted.mean(x = ar, w = weight),
    .groups = "keep"
  ) %>%
  ungroup() %>%
  # redact if <=5 events
  mutate(across(c("events", "ar_crude", "ar_weighted"), 
                ~if_else(events <= 5, NA_real_, as.numeric(.x)))) %>%
  select(-events)

readr::write_csv(
  data_ar_arm_weight,
  here::here("output", "tte", "tables", "ar_arm_weight.csv")
)

### 3. age- and sex-standardised AR within cohort
# AR grouped by age, sex
data_counts <- data_tte %>%
  left_join(data_all, by = "patient_id") %>% 
  group_by(cohort, outcome, k, age, sex) %>%
  summarise(
    n = n(),
    events = sum(status),
    .groups = "keep"
  ) %>%
  ungroup() %>%
  mutate(ar = events/n) 

# within each outcome and comparison period, 
# get proportion of general population in each sex:age_band
arm_weights <- data_counts %>%
  filter(cohort == "noncancer") %>%
  group_by(outcome, k) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  mutate(prop = n/total) %>%
  select(outcome, k, age, sex, weight = prop)

# calculate AR in both arms, 
# weighted by sex:age_band distribution in noncancer pop
data_ar <- data_counts %>%
  left_join(
    arm_weights, by = c("outcome", "k", "age", "sex")
  ) %>%
  group_by(cohort, outcome, k) %>%
  summarise(
    events = sum(events),
    ar_crude = sum(events)/sum(n),
    ar_weighted = weighted.mean(x = ar, w = weight),
    .groups = "keep"
  ) %>%
  ungroup() %>%
  # redact if <=5 events
  mutate(across(c("events", "ar_crude", "ar_weighted"), 
                ~if_else(events <= 5, NA_real_, as.numeric(.x)))) %>%
  select(-events)

readr::write_csv(
  data_ar,
  here::here("output", "tte", "tables", "ar_agesexonly.csv")
)
