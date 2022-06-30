library(tidyverse)
library(glue)

################################################################################
# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))


# eligibility_count_all <- readr::read_csv(
#   here::here(release_folder, "eligibility_count_all.csv")
# )
# 
# eligibility_count_all %>%
#   mutate(across(c(n,n_removed), ~scales::comma(.x, accuracy = 1))) %>%
#   print(n=Inf)
# 
# #TODO
# # sort out boxes E and pink for both vax and unvax

################################################################################
table1_main_TRUE_REDACTED <- readr::read_csv(
  "release20220622/table1_main_TRUE_REDACTED.csv"
  )
cols <- names(table1_main_TRUE_REDACTED[,3:8])
counts <- as.numeric(str_remove_all((unlist(table1_main_TRUE_REDACTED[1,3:8])), ","))
names(counts) <- cols
# cancer cohort:
scales::comma(sum(counts[str_detect(names(counts), "_1")]))
# general cohort
scales::comma(sum(counts[str_detect(names(counts), "_2")]))


################################################################################
event_counts <- readr::read_csv(
  here::here(release_folder, "event_counts.csv")
)

event_counts %>%
  filter(outcome == "coviddeath", k==1) %>%
  group_by(subgroup) %>%
  summarise(total = scales::comma(sum(n)))



################################################################################
estimates_all <- readr::read_csv(
  here::here(release_folder, "estimates_all.csv")
) %>%
  mutate(across(c(estimate, conf.low, conf.high), exp)) 

estimates_all %>%
  filter(
    variable == "k",
    subgroup %in% 1:2,
    !reference_row,
    model=="unadjusted",
    outcome %in% c("postest", "covidadmitted", "coviddeath"),
    comparison=="both",
    prior
    ) %>%
  group_by(subgroup, outcome) %>%
  mutate(min_k = min(period), max_k = max(period)) %>%
  ungroup() %>%
  filter(period==min_k|period==max_k) %>%
  mutate(across(c(estimate, conf.low, conf.high), ~round(100*(1-.x),1)))  %>%
  transmute(
    subgroup, outcome, period, 
    value = glue("{estimate}% ({conf.low}-{conf.high})")
  )


# percentage point differences in 3-6 o 23-26
ppoint <- estimates_all %>%
  filter(
    variable == "k",
    subgroup %in% 3:4,
    !reference_row,
    model=="unadjusted",
    outcome %in% c("postest", "covidadmitted", "coviddeath"),
    comparison=="both",
    prior
    ) %>%
  mutate(across(c(estimate, conf.low, conf.high), ~round(100*(1-.x),1)))  

ppoint %>%
  select(subgroup, outcome, period, estimate) %>%
  pivot_wider(
    names_from = period,
    values_from = estimate
  ) %>%
  mutate(
    `1` - `6`,
    `2` - `6`
    )



metareg_results_rhr <- readr::read_rds(
  here::here(release_folder, "metareg_results_rhr.rds")
)

metareg_results_rhr %>% 
  filter(
    subgroup %in% subgroups[3:4],
    comparison=="both", 
    outcome%in%c("postest","covidadmitted", "coviddeath"),
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
