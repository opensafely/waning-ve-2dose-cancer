###############################################################################
library(tidyverse)
library(glue)

###############################################################################
# tabulate event counts for each of the analyses

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

event_counts <- data_tte %>%
  group_by(analysis, subgroup, arm, outcome, k) %>%
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

# save for release
readr::write_csv(
  event_counts,
  here::here("output", "tte", "tables", "event_counts.csv")
)