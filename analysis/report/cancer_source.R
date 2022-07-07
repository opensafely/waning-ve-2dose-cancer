################################################################################
# summarise sources of cancer data

################################################################################
library(tidyverse)

fs::dir_create(here::here("output", "report"))

# read data for ever covariates
data_cancer <- arrow::read_feather(
  file = here::here("output", "input_covs.feather")) %>%
  # all cancer defined on or after 2018-01-01
  transmute(
    patient_id,
    cancer_icd10 = pmax(
      !is.na(cancer_solid_icd10_date),
      !is.na(cancer_haem_icd10_date),
      !is.na(cancer_unspec_icd10_date)
      ),
    cancer_snomed = pmax(
      !is.na(cancer_solid_snomed_date),
      !is.na(cancer_haem_snomed_date)
    )
  ) %>%
  mutate(across(starts_with("cancer"), as.logical))

capture.output(
  data_cancer %>%
    group_by(cancer_snomed, cancer_icd10) %>%
    count() %>%
    ungroup() %>%
    mutate(percent = round(100*n/sum(n), 3)) %>%
    kableExtra::kable(
      format = "pipe",
      caption = "Cross-tab of cancer identified from SNOMED and ICD10 codes"
    ),
  file = here::here("output", "report", "cancer_source.txt"),
  append = FALSE
)
