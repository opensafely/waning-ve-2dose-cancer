################################################################################

# This script:
# - reads, processes the extracted data, saves the following:
# - data_processed.rds = processed data to be used for applying eligibility criteria
# - data_processed_0_tabulate.txt = summary of extracted data for checking
# - data_*_vax_dates.rds = long and wide vaccine data

################################################################################

## setup
library(tidyverse)
library(lubridate)
library(glue)

## source functions
source(here::here("analysis", "functions", "data_process_functions.R"))
source(here::here("analysis", "functions", "data_properties.R"))

## create folders for outputs
fs::dir_create(here::here("output", "data"))
fs::dir_create(here::here("output", "tables"))

## import study_parameters
study_parameters <- readr::read_rds(
  here::here("analysis", "lib", "study_parameters.rds"))

## regions
regions <- readr::read_csv(
  here::here("analysis", "lib", "regions.csv")
)

################################################################################
# initial pre-processing
cat("#### extract data ####\n")

data_extract <- 
  arrow::read_feather(file = here::here("output", "input_vax.feather")) %>%
  # because date types are not returned consistently by cohort extractor
  mutate(across(c(contains("_date")), 
                ~ floor_date(
                  as.Date(., format="%Y-%m-%d"),
                  unit = "days"))) %>%
  mutate(across(imd, ~as.integer(as.character(.x))))

cat("#### process extracted data ####\n")
data_processed_0 <- data_extract %>%
  # derive ethnicity variable
  mutate(
    # Region
    region = factor(region, levels = regions$region),
    # Ethnicity
    ethnicity = if_else(is.na(ethnicity_6), ethnicity_6_sus, ethnicity_6),
    ethnicity = fct_case_when(
      ethnicity == "1" ~ "White",
      ethnicity == "4" ~ "Black",
      ethnicity == "3" ~ "South Asian",
      ethnicity == "2" ~ "Mixed",
      ethnicity == "5" ~ "Other",
      TRUE ~ NA_character_
    ),
    # IMD quintile
    imd = fct_case_when(
      imd < 1 | is.na(imd) ~ NA_character_,
      imd < 32844*1/5 ~ "1 most deprived",
      imd < 32844*2/5 ~ "2",
      imd < 32844*3/5 ~ "3",
      imd < 32844*4/5 ~ "4",
      TRUE ~ "5 least deprived"
    ),
    # Sex
    sex = fct_case_when(
      sex == "F" ~ "Female",
      sex == "M" ~ "Male",
      TRUE ~ NA_character_
    ),
    #Subgroup
    subgroup = fct_case_when(
      jcvi_group %in% c("04b", "06") & age_1 < 65 ~ "18-64 years and clinically vulnerable",
      jcvi_group %in% c("11", "12") ~ "18-39 years",
      jcvi_group %in% c("07", "08", "09", "10") ~ "40-64 years",
      jcvi_group %in% c("02", "03", "04a", "04b", "05") ~ "65+ years",
      TRUE ~ NA_character_
    )
    
  ) %>%
  select(-ethnicity_6, -ethnicity_6_sus) %>%
  droplevels()

################################################################################
cat("#### properties of data_processed_0 ####\n")
# for checking for errors
data_properties(
  data = data_processed_0,
  path = file.path("output", "tables")
)

cat("## check subgroups as desired ##\n")
data_processed_0 %>%
  group_by(subgroup, jcvi_group) %>%
  count() %>%
  ungroup()

################################################################################
# process vaccine data
data_vax <- local({
  
  data_vax_pfizer <- data_processed_0 %>%
    select(patient_id, matches("covid\\_vax\\_pfizer\\_\\d+\\_date")) %>%
    pivot_longer(
      cols = -patient_id,
      names_to = c(NA, "vax_pfizer_index"),
      names_pattern = "^(.*)_(\\d+)_date",
      values_to = "date",
      values_drop_na = TRUE
    ) %>%
    arrange(patient_id, date)
  
  data_vax_az <- data_processed_0 %>%
    select(patient_id, matches("covid\\_vax\\_az\\_\\d+\\_date")) %>%
    pivot_longer(
      cols = -patient_id,
      names_to = c(NA, "vax_az_index"),
      names_pattern = "^(.*)_(\\d+)_date",
      values_to = "date",
      values_drop_na = TRUE
    ) %>%
    arrange(patient_id, date)
  
  data_vax_moderna <- data_processed_0 %>%
    select(patient_id, matches("covid\\_vax\\_moderna\\_\\d+\\_date")) %>%
    pivot_longer(
      cols = -patient_id,
      names_to = c(NA, "vax_moderna_index"),
      names_pattern = "^(.*)_(\\d+)_date",
      values_to = "date",
      values_drop_na = TRUE
    ) %>%
    arrange(patient_id, date)
  
  data_vax_disease <- data_processed_0 %>%
    select(patient_id, matches("covid\\_vax\\_disease\\_\\d+\\_date")) %>%
    pivot_longer(
      cols = -patient_id,
      names_to = c(NA, "vax_disease_index"),
      names_pattern = "^(.*)_(\\d+)_date",
      values_to = "date",
      values_drop_na = TRUE
    ) %>%
    arrange(patient_id, date)
  
  
  data_vax <- data_processed_0 %>% # to get the unvaccinated
    # filter(if_all(starts_with("covid_vax"), ~ is.na(.))) %>%
    filter_at(vars(starts_with("covid_vax")), all_vars(is.na(.))) %>%
    select(patient_id) %>% 
    full_join(
      data_vax_pfizer %>%
        full_join(data_vax_az, by=c("patient_id", "date")) %>%
        full_join(data_vax_moderna, by=c("patient_id", "date")) %>%
        full_join(data_vax_disease, by=c("patient_id", "date")),
      by = "patient_id"
    ) %>%
    mutate(
      brand = fct_case_when(
        (!is.na(vax_az_index)) & is.na(vax_pfizer_index) & is.na(vax_moderna_index) ~ "az",
        is.na(vax_az_index) & (!is.na(vax_pfizer_index)) & is.na(vax_moderna_index) ~ "pfizer",
        is.na(vax_az_index) & is.na(vax_pfizer_index) & (!is.na(vax_moderna_index)) ~ "moderna",
        (!is.na(vax_az_index)) + (!is.na(vax_pfizer_index)) + (!is.na(vax_moderna_index)) > 1 ~ "duplicate",
        !is.na(vax_disease_index) ~ "unknown",
        TRUE ~ NA_character_
      )
    ) %>%
    arrange(patient_id, date) %>%
    group_by(patient_id) %>%
    mutate(
      vax_index=row_number()
    ) %>%
    ungroup() %>%
    droplevels()
  
  data_vax
  
})

data_vax_wide <- data_vax %>%
  pivot_wider(
    id_cols= patient_id,
    names_from = c("vax_index"),
    values_from = c("date", "brand"),
    names_glue = "covid_vax_{vax_index}_{.value}"
  )

readr::write_rds(
  data_vax_wide,
  here::here("output", "data", "data_wide_vax_dates.rds"), 
  compress="gz")

################################################################################
# save dataset of covariates 
# (remove variables that are saved elsewhere)

data_processed <- data_processed_0 %>%
  select(
    # remove vaccine variables
    -contains("_vax_"),
  ) 

readr::write_rds(
  data_processed,
  here::here("output", "data", "data_processed.rds"), 
  compress="gz")
