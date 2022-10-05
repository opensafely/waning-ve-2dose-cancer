################################################################################
# This script:
# - combines estimates from all models into csv for release

################################################################################
library(tidyverse)
library(RColorBrewer)
library(glue)

################################################################################
fs::dir_create(here::here("output", "release_objects"))

################################################################################
# read study parameters
study_parameters <- readr::read_rds(
  here::here("analysis", "lib", "study_parameters.rds"))
K <- study_parameters$K

# read outcomes
outcomes <- readr::read_rds(
  here::here("analysis", "lib", "outcomes.rds")
)

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)

# define comparisons
comparisons <- readr::read_rds(
  here::here("analysis", "lib", "comparisons.rds"))

# redaction functions
source(here::here("analysis", "functions", "redaction_functions.R"))

################################################################################
# model estimates
list_modelcox_tidy <- function(
  c = comparisons,
  s = subgroup_labels, 
  p = TRUE
) {
  unlist(
    lapply(
      c,
      function(x)
        lapply(
          s,
          function(y)
            list.files(path = here::here("output", "models_cox", "data"),
                       pattern = glue("modelcox_tidy_{x}_{y}_{p}_.+.rds"),
                       all.files = FALSE,
                       full.names = FALSE, recursive = FALSE,
                       ignore.case = FALSE, include.dirs = FALSE)
        )
      
    ) 
  )
}


all_files <- c(
  list_modelcox_tidy(),
  list_modelcox_tidy(s=1:2,p=FALSE)
)

mtimes <- sapply(
  all_files,
  function(x)
    format(file.info(here::here("output", "models_cox", "data", x))$mtime, "%m/%d/%y %H:%M")
)

cat("Check when files were last modified:")
tibble(
  file = names(mtimes),
  modified = unname(mtimes)
) %>%
  arrange(modified) %>%
  print(n=Inf)

# read files
model_tidy_list <- lapply(
  all_files,
  function(filename) {
    filename_split <- unlist(str_split(str_remove(filename, ".rds"), "_"))
    readr::read_rds(
      here::here("output", "models_cox", "data", filename)
    ) %>%
      mutate(
        comparison = filename_split[3],
        subgroup = filename_split[4],
        prior = as.logical(filename_split[5]),
        outcome = filename_split[6],
        period = filename_split[7]
        )
  }
)

model_tidy_tibble <- bind_rows(
  model_tidy_list[sapply(model_tidy_list, function(x) is_tibble(x))]
) %>%
  # mutate(across(label, 
  #               ~ if_else(
  #                 variable == "k" & label != "0",
  #                 period,
  #                 .x
  #               ))) %>%
  mutate(across(c(estimate, conf.low, conf.high), round, 5)) %>%
  mutate(across(model, 
                factor, levels = 1:3, labels = c("unadjusted", "part_adjusted", "max_adjusted"))) %>%
  # calculate the total number of observations per model
  mutate(n_obs_model = if_else(variable == "k", n_obs, NA_real_)) %>%
  group_by(subgroup, comparison, outcome, model, period) %>%
  mutate(across(n_obs_model, sum, na.rm=TRUE)) %>%
  ungroup() %>%
  # round n_obs_model up to nearest 7
  mutate(across(c(n_obs_model, exposure), ~ceiling_any(.x, to=7))) %>%
  select(subgroup, comparison, prior, outcome, model, period, variable, label, reference_row,
         n_obs_model, estimate, conf.low, conf.high) 

# check size of smallest model
print(min(model_tidy_tibble$n_obs_model, na.rm=TRUE))

readr::write_csv(
  model_tidy_tibble,
  here::here("output", "release_objects", "estimates_all.csv"))
