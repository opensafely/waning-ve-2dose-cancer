################################################################################

# This script:
# creates time-to-event data for all outcomes

################################################################################

library(tidyverse)
library(glue)

################################################################################
study_parameters <- readr::read_rds(
  here::here("analysis", "lib", "study_parameters.rds"))
K <- study_parameters$K

# read outcomes
outcomes <- readr::read_rds(
  here::here("analysis", "lib", "outcomes.rds"))
outcomes_death <- outcomes[str_detect(outcomes, "death")]

################################################################################
# read data
data_all <- readr::read_rds(
  here::here("output", "data", "data_all.rds")) 

################################################################################
# redaction functions
source(here::here("analysis", "functions", "redaction_functions.R"))

################################################################################
# output directories
fs::dir_create(here::here("output", "tte", "data"))

################################################################################

data <- data_all %>%
  # reshape such that each line corresponds to a comparison period
  pivot_longer(
    cols = matches("\\w+_\\d_date"),
    names_to = c(".value", "k"),
    names_pattern = "(.*)_(.)_date"
  ) %>%
  rename_with(~str_c(.x, "_date"), .cols = all_of(c("start", "end", "anytest"))) %>%
  mutate(across(k, as.integer)) %>%
  # keep only odd unvax for odd k, equiv. for even
  # see original waning protocol for rational for splitting the unvaccinated group
  filter(
    is.na(split) |
      ((k %% 2) == 0 & split == "even") |
      ((k %% 2) != 0 & split == "odd")
  ) %>%
  select(patient_id, k, arm, ends_with("date"), ends_with("subgroup"), starts_with("exclude"))

################################################################################
# generates and saves data_tte and tabulates event counts 
# returns tables of events
derive_data_tte <- function(
  .data, 
  outcome
  ) {
  
  # remove periods for which outcome has occurred before the start date
  # e.g. 
  # if outcome is postest and postest occurs in period 3, remove periods 4,5,6;
  # if outcome is anytest and anytest occurs in period 3, keep periods 4,5,6
  # BUT if outcome is anytest and postest occurs in period 3, remove periods 4,5,6
  if (outcome == "anytest") {
    outcome_exclude <- "postest"
  } else {
    outcome_exclude <- outcome
  }
  
  # function to be applied in dplyr::filter
  occurs_after_start_date <- function(cov_date, index_date) {
    is.na(cov_date) | index_date < cov_date
  }
  
  data_tte <- .data %>%
    # exclude if subsequent_vax, death, dereg or outcome_exclude occurred before start of period
    filter_at(
      vars(str_c(unique(c("subsequent_vax", "dereg", "coviddeath", "noncoviddeath", outcome_exclude)), "_date")),
      all_vars(occurs_after_start_date(cov_date = ., index_date = start_date))
    ) %>%
    # only keep periods for which start_date < study end_date
    filter(
      start_date < as.Date(study_parameters$end_date) 
    ) %>%
    # if end_date > study_parameters$end_date, replace with study_parameters$end_date
    mutate(across(end_date,
                  ~ if_else(as.Date(study_parameters$end_date) < .x,
                            as.Date(study_parameters$end_date),
                            .x))) %>%
    # only keep dates for censoring and outcome variables between start_date and end_date
    mutate(across(all_of(str_c(unique(c("dereg", outcomes)), "_date")),
                  ~ if_else(
                    !is.na(.x) & (start_date < .x) & (.x <= end_date),
                    .x,
                    as.Date(NA_character_)
                  ))) %>%
    # new time-scale: time since earliest start_fu_date in data
    mutate(across(ends_with("_date"),
                  ~ as.integer(.x - min(start_date)))) %>%
    rename_at(vars(ends_with("_date")),
              ~ str_remove(.x, "_date")) %>%
    mutate(
      # censor follow-up time at first of the following:
      tte = pmin(!! sym(outcome), dereg, coviddeath, noncoviddeath, end, na.rm = TRUE),
      status = if_else(
        !is.na(!! sym(outcome)) & !! sym(outcome) == tte,
        TRUE,
        FALSE
      )) %>%
    select(patient_id, arm, k, tstart = start, tstop = tte, status,
           ends_with("subgroup"), starts_with("exclude")) %>%
    arrange(patient_id, k) %>%
    mutate(tmp_outcome = outcome) %>%
    # update prior_infection_subgroup for specific outcome variables
    # exclude_postest = postest in [start_date - 27, start_date]
    # exclude_covidadmitted = covidadmitted in [start_date - 27, start_date]
    mutate(across(prior_infection_subgroup,
                  ~ case_when(
                    !.x & 
                      !(tmp_outcome %in% c("postest", "anytest") & exclude_postest) & 
                      !(tmp_outcome == "covidadmitted" & exclude_covidadmitted) ~ FALSE,
                    TRUE ~ TRUE
                  ))) %>%
    select(-tmp_outcome, -starts_with("exclude_")) %>%
    mutate(across(arm, factor, levels = c("BNT162b2", "ChAdOx1", "unvax"))) %>%
    mutate(across(k, factor, levels = 1:K)) 
  
  # checks
  stopifnot("tstart should be  >= 0 in data_tte" = data_tte$tstart>=0)
  stopifnot("tstop - tstart should be strictly > 0 in data_tte" = data_tte$tstop - data_tte$tstart > 0)
  
  # save data_tte
  readr::write_rds(
    data_tte,
    here::here("output", "tte", "data", glue("data_tte_{outcome}.rds")),
    compress = "gz")
  
}

################################################################################
# apply derive_data_tte for all outcomes

for (x in outcomes) {
  data %>% derive_data_tte(outcome = x)
}


