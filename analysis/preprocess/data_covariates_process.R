################################################################################
# process comparisons data
library(tidyverse)
library(lubridate)
library(glue)

################################################################################
## import study_parameters
study_parameters <- readr::read_rds(
  here::here("analysis", "lib", "study_parameters.rds"))
K <- study_parameters$K

# import variable names
model_varlist <- readr::read_rds(
  here::here("analysis", "lib", "model_varlist.rds")
)

# read outcomes
outcomes <- readr::read_rds(
  here::here("analysis", "lib", "outcomes.rds"))

################################################################################
# individuals eligible based on box c, d & e criteria 
# arm and split info
data_arm <- bind_rows(
  readr::read_rds(
    here::here("output", "data", "data_eligible_e_vax.rds")) %>%
    rename(arm=brand),
  readr::read_rds(
    here::here("output", "data", "data_eligible_e_unvax.rds")) %>%
    mutate(arm = "unvax")
)  %>%
  select(patient_id, arm, split)

# vars from data_processed
data_processed <- readr::read_rds(
  here::here("output", "data", "data_processed.rds")) %>%
  select(patient_id, #subgroup,
         jcvi_group, elig_date, region, ethnicity,
         dereg_date, death_date,
         any_of(unname(model_varlist$demographic)))

# vax data
data_wide_vax_dates <- readRDS(
  here::here("output", "data", "data_wide_vax_dates.rds")) %>%
  select(patient_id, covid_vax_1_date, covid_vax_3_date)

# read data for ever covariates
data_covariates <- arrow::read_feather(
  file = here::here("output", "input_covs.feather")) 

################################################################################
data_all <- data_arm %>%
  # join to covariates data
  left_join(
    data_covariates %>%
      select(patient_id, 
             matches(c("start_\\d_date", "end_\\d_date")),
             starts_with("cancer"),
             starts_with(c(unname(outcomes), "primary_care_covid")),
             any_of(unname(unlist(model_varlist)))) %>%
      mutate(across(contains("_date"), 
                    ~ floor_date(
                      as.Date(.x, format="%Y-%m-%d"),
                      unit = "days"))),
    by = "patient_id") %>%
  # join to data_processed
  left_join(
    data_processed, by = "patient_id"
  ) %>%
  # join to vaccines
  left_join(
    data_wide_vax_dates, 
    by = "patient_id"
  ) %>%
  # derive remaining covariates
  mutate(
    
    pregnancy = pregnancy & (sex == "Female") & (age < 50),
    
    multimorb =
      # as.integer(bmi %in% "Obese III (40+)") +
      as.integer(chd)  +
      as.integer(diabetes) +
      as.integer(cld) +
      as.integer(ckd) +
      as.integer(crd) +
      as.integer(immunosuppressed) +
      as.integer(cns),
    
    multimorb = cut(
      multimorb,
      breaks = c(0, 1, 2, Inf),
      labels=c("0", "1", "2+"),
      right=FALSE),
    
    noncoviddeath_date = if_else(
      !is.na(death_date) & is.na(coviddeath_date),
      death_date,
      as.Date(NA_character_)
      )
    
  ) %>%
  mutate(across(test_hist_n,
                ~ factor(case_when(
                  is.na(.x) ~ NA_character_,
                  .x < 1 ~ "0",
                  .x < 2 ~ "1",
                  .x < 3 ~ "2",
                  TRUE ~ "3+"
                )))) %>%
  mutate(subsequent_vax_date = if_else(
    arm %in% "unvax",
    covid_vax_1_date,
    covid_vax_3_date)) %>%
  select(-covid_vax_1_date, -covid_vax_3_date) %>%
  mutate(
    
    # Create flags for each subgroup analysis:
    # positive test, primary care, or hospitalisation on or before start_1_date - 28 days
    prior_infection_subgroup = if_else(
      !is.na(postest_0_date) | !is.na(primary_care_covid_case_0_date) | !is.na(covidadmitted_0_date),
      TRUE, FALSE
    ),
    
    # those to exclude when positive test is outcome
    exclude_postest = if_else(
      !is.na(postest_1_date) | !is.na(primary_care_covid_case_1_date) | !is.na(covidadmitted_1_date),
      TRUE, FALSE
    ),
    
    # those to exclude when hospitalisation outcome
    exclude_covidadmitted = if_else(
      !is.na(covidadmitted_1_date),
      TRUE, FALSE
    ),
    
    # cancer subgroups
    cancer_subgroup = case_when(
      !is.na(cancer_haem_icd10_date) | !is.na(cancer_haem_snomed_date) ~ "cancer_haem",
      !is.na(cancer_solid_icd10_date) | !is.na(cancer_solid_snomed_date) ~ "cancer_solid",
      !is.na(cancer_unspec_icd10_date) ~ "cancer_unspec",
      TRUE ~ "noncancer"
    ),
    
    # age subgroups
    age_subgroup = factor(
      if_else(age < 70, "18-69", "70+"),
      levels = c("18-69", "70+")
    )
  ) %>%
  rename(
    postest_date = postest_2_date,
    covidadmitted_date = covidadmitted_2_date
    ) %>%
  select(
    -matches(c("cancer_.+_date", "postest_\\d_date", "covidadmitted_\\d_date", "primary_care_covid_case_\\d_date")),
    -any_of(unname(model_varlist$multimorb))
  )

readr::write_rds(
  data_all,
  here::here("output", "data", "data_all.rds"),
  compress = "gz"
)
