library(tidyverse)

if (!exists("release_folder")) release_folder <- here::here("output", "release_objects")

# read data
estimates_all <- readr::read_csv(file.path(release_folder, "estimates_all.csv"))

# read outcomes
outcomes <- readr::read_rds(
  here::here("analysis", "lib", "outcomes.rds")
)

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))

# process
data_metareg_0 <- estimates_all %>%
  filter(
    outcome %in% outcomes,
    !reference_row,
    variable %in% "k"
    ) %>%
  mutate(
    subgroup = as.integer(str_extract(subgroup, "\\d"))
    ) %>%
  select(subgroup, comparison, outcome, model, k = label, estimate, conf.low, conf.high) 

# expand to all combinations
data_metareg <- data_metareg_0 %>%
  expand(subgroup, comparison, outcome, model, k) %>%
  left_join(data_metareg_0) %>%
  # following lines because in a few instances seloghr=0, and this is breaking the metareg code 
  mutate(seloghr = (conf.high - conf.low)/(2*1.96)) %>%
  mutate(across(c(estimate, conf.low, conf.high),
                ~if_else(near(seloghr, 0), NA_real_, .x))) %>%
  select(-seloghr)

# save data
readr::write_csv(
  data_metareg,
  here::here(release_folder, "data_metareg.csv"))
