library(tidyverse)

release_folder <- here::here("release20220622")

metareg_results_rhr <- readr::read_rds(
  here::here(release_folder, "metareg_results_rhr.rds")
)

################################################################################
# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))

################################################################################


for (m in 1:3) {
  
  ftab <- metareg_results_rhr %>%
    filter(
      model==m,
      prior,
      outcome != "anytest"
    ) %>%
    mutate(across(starts_with("rhr"), ~round(.x,2))) %>%
    transmute(
      outcome, subgroup, comparison,
      value = glue::glue("{format(rhr, nsmall=2)} ({format(rhr_lower,nsmall=2)}, {format(rhr_higher, nsmall=2)})")
    ) %>%
    pivot_wider(
      names_from = outcome,
      values_from = value
    ) %>%
    select(
      Subgroup = subgroup, Brand = comparison, 
      "Positive SARS-CoV-2 test" = postest,
      "COVID-19 hospitalisation" = covidadmitted,
      "COVID-19 death" = coviddeath,
      "non-COVID-19 death" = noncoviddeath,
    ) %>%
    mutate(across(Subgroup,
                  factor, 
                  levels = subgroups[c(1:2, 5:8, 3:4)],
                  labels = c(
                    "Cancer", "Non-cancer", 
                    "Cancer, 18-69 years", "Cancer, 70+", 
                    "Non-cancer, 18-69 years", "Non-cancer, 70+",
                    "Haematological malignancy", "Solid organ malignancy"
                  ))) %>%
    mutate(across(Brand,
                  factor,
                  levels = c("BNT162b2", "ChAdOx1", "both"),
                  labels = c("BNT162b2", "ChAdOx1", "Combined"))) %>%
    arrange(Subgroup, Brand) %>%
    flextable::flextable() %>%
    flextable::merge_v(j=1, part = "body") %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::theme_booktabs() 
  
  doc <- officer::read_docx() %>%
    flextable::body_add_flextable(value = ftab, split = FALSE) %>%
    print(target = here::here(release_folder, glue::glue("rhr_table_{m}.docx")))
  
}

