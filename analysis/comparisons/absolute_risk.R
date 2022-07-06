###############################################################################
library(tidyverse)
library(glue)

###############################################################################
# age- and sex- standardised absolute risk

###############################################################################
# read outcomes
outcomes <- readr::read_rds(
  here::here("analysis", "lib", "outcomes.rds"))

################################################################################
# redaction functions
source(here::here("analysis", "functions", "redaction_functions.R"))

################################################################################
# create output directory
fs::dir_create(here::here("output", "tte", "tables"))

###############################################################################
# read tte data
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
    subgroup = if_else(
      cancer_subgroup %in% "noncancer", 
      cancer_subgroup, 
      "cancer"
      )
    ) %>%
  mutate(across(arm, ~if_else(.x %in% "unvax", as.character(.x), "vax"))) %>%
  select(patient_id, subgroup, arm, k, outcome, status)

# read age and sex data
age_bands <- c(-Inf, seq(25,90,5), Inf)
lower_bands <- as.character(age_bands[-length(age_bands)])
lower_bands[1] <- "18"
upper_bands <- as.character(age_bands[-1] - 1)
upper_bands[length(upper_bands)] <- "120"
age_band_labs <- str_c(lower_bands, upper_bands, sep = "-")

data_all <- readr::read_rds(
  here::here("output", "data", "data_all.rds")) %>%
  select(patient_id, age, sex) %>%
  mutate(across(age, ~cut(.x, breaks = age_bands, labels = age_band_labs)))

###############################################################################
# risk grouped by arm
data_counts_arm_group <- data_tte %>%
  left_join(data_all, by = "patient_id") %>% 
  group_by(subgroup, arm, outcome, k, age, sex) %>%
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
  filter(subgroup == "noncancer") %>%
  group_by(outcome, arm, k) %>%
  mutate(prop = n/sum(n)) %>%
  # mutate(total = sum(n)) %>%
  ungroup() %>%
  select(arm, outcome, k, age, sex, weight = prop)

# calculate ar in both groups, weighted by sex:age_band distribution in noncancer pop
data_ar_arm_group <- data_counts_arm_group %>%
  left_join(
    weights_arm_group, by = c("arm", "outcome", "k", "age", "sex")
  ) %>%
  group_by(subgroup, arm, outcome, k) %>%
  summarise(
    ar_crude = sum(events)/sum(n),
    ar_weighted = weighted.mean(x = ar, w = weight),
    .groups = "keep"
  ) %>%
  ungroup() 

readr::write_csv(
  data_ar_arm_group,
  here::here("output", "tte", "tables", "ar_arm_group.csv")
  )

###
# don't group by arm, but include in weights

# within each outcome and comparison period, 
# get proportion of general population in each arm:sex:age_band
weights_arm_weight <- data_counts_arm_group %>%
  filter(subgroup == "noncancer") %>%
  group_by(outcome, k) %>%
  mutate(prop = n/sum(n)) %>%
  # mutate(total = sum(n)) %>%
  ungroup() %>%
  select(arm, outcome, k, age, sex, weight = prop)

# calculate ar in both groups, weighted by sex:age_band distribution in noncancer pop
data_ar_arm_weight <- data_counts_arm_group %>%
  left_join(
    weights_arm_weight, by = c("arm", "outcome", "k", "age", "sex")
  ) %>%
  group_by(subgroup, outcome, k) %>%
  summarise(
    ar_crude = sum(events)/sum(n),
    ar_weighted = weighted.mean(x = ar, w = weight),
    .groups = "keep"
  ) %>%
  ungroup() 

readr::write_csv(
  data_ar_arm_weight,
  here::here("output", "tte", "tables", "ar_arm_weight.csv")
)

###
# don't group or weight by arm
data_counts <- data_tte %>%
  left_join(data_all, by = "patient_id") %>% 
  group_by(subgroup, outcome, k, age, sex) %>%
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
  filter(subgroup == "noncancer") %>%
  group_by(outcome, k) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  mutate(prop = n/total) %>%
  select(outcome, k, age, sex, weight = prop)

# calculate ar in both groups, weighted by sex:age_band distribution in noncancer pop
data_ar <- data_counts %>%
  left_join(
    arm_weights, by = c("outcome", "k", "age", "sex")
  ) %>%
  group_by(subgroup, outcome, k) %>%
  summarise(
    ar_crude = sum(events)/sum(n),
    ar_weighted = weighted.mean(x = ar, w = weight),
    .groups = "keep"
  ) %>%
  ungroup() 

readr::write_csv(
  data_ar,
  here::here("output", "tte", "tables", "ar_agesexonly.csv")
)
