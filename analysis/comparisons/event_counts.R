###############################################################################
library(tidyverse)
library(glue)

###############################################################################
# tabulate event counts for each of the subgroups

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
  mutate(persondays = tstop - tstart) 
  
###############################################################################
# event counts function

event_count_function <- function(.data, include_prior_covid = TRUE) {
  
  if (include_prior_covid) {
    .data <- .data
  } else {
    .data <- .data %>% 
      filter(!prior_infection_subgroup) %>%
      mutate(across(subgroup, ~ str_c(.x, "_no_prior_covid")))
  }
  
  .data %>%
    group_by(subgroup, arm, outcome, k) %>%
    summarise(
      # personyears rounded to nearest whole year
      personyears = round(sum(persondays)/365.25,0),
      # individuals rounded up to the nearest 7
      n = ceiling_any(n(), to = 7),
      # events rounded up to the nearest 7
      events = ceiling_any(sum(status), to = 7), 
      .groups = "keep"
    ) %>%
    ungroup() 
  
}

# main analysis
event_counts_main <- bind_rows(lapply(
  c(TRUE, FALSE),
  function(x) 
    data_tte %>%
    # derive cancer / nocancer subgroup
    mutate(
      subgroup = if_else(cancer_subgroup == "noncancer", "noncancer", "cancer")
    ) %>%
    event_count_function(include_prior_covid = x)
))

# cancer type analysis
event_counts_type <- bind_rows(lapply(
  c(TRUE, FALSE),
  function(x) 
    data_tte %>%
    # keep only known haem and solid organ cancer subgroups
    filter(cancer_subgroup %in% c("cancer_haem", "cancer_solid")) %>%
    mutate(subgroup = str_remove(cancer_subgroup, "cancer_")) %>%
    event_count_function(include_prior_covid = x)
))

# age analysis
event_counts_age <- bind_rows(lapply(
  c(TRUE, FALSE),
  function(x) 
    data_tte %>%
    # derive (cancer / nocancer) * (18-69 / 70+) subgroups
    mutate(
      subgroup = if_else(cancer_subgroup == "noncancer", "noncancer", "cancer")
    ) %>%
    mutate(across(subgroup, ~ str_c(subgroup, "_", age_subgroup))) %>%
    event_count_function(include_prior_covid = x)
))

# bind all event type tables
event_counts <- bind_rows(event_counts_main, event_counts_type, event_counts_age)
  
# save for release
readr::write_csv(
  event_counts,
  here::here("output", "tte", "tables", "event_counts.csv")
)