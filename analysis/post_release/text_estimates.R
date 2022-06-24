

metareg_results_rhr <- readr::read_rds(
  here::here(release_folder, "metareg_results_rhr.rds")
)

metareg_results_rhr %>% 
  filter(
    subgroup%in% subgroups[1:2],
    comparison=="both", 
    outcome=="postest",
    prior,
    model==1
    ) %>%
  mutate(across(starts_with("rhr"), ~round(.x,2))) %>%
  mutate(value = glue("{rhr} [{rhr_lower}-{rhr_higher}]"))

metareg_results_rhr %>% 
  filter(
    subgroup%in% subgroups[1:2],
    comparison=="both", 
    outcome%in%c("covidadmitted", "coviddeath"),
    prior,
    model==1
  ) %>%
  mutate(across(starts_with("rhr"), ~round(.x,2))) %>%
  mutate(value = glue("{rhr} [{rhr_lower}-{rhr_higher}]"))
