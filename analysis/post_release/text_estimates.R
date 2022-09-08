library(tidyverse)
library(glue)

release_folder <- here::here("release20220622")

################################################################################
# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))

################################################################################
# Study population results section

table1_main_TRUE_REDACTED <- readr::read_csv(
  "release20220622/table1_main_TRUE_REDACTED.csv"
  )
cols <- names(table1_main_TRUE_REDACTED[,3:8])
counts <- as.numeric(str_remove_all((unlist(table1_main_TRUE_REDACTED[1,3:8])), ","))
names(counts) <- cols

## number eligibile for each group
# BNT162b2
sum(counts[c(1,4)])
# ChAdOx1
sum(counts[c(2,5)])
# unvax
sum(counts[c(3,6)])

## percent of each group in the cancer cohort
# BNT162b2
100*counts[1]/sum(counts[c(1,4)])
# ChAdOx1
100*counts[2]/sum(counts[c(2,5)])
# unvax
100*counts[3]/sum(counts[c(3,6)])

## total in each cohort
# cancer cohort:
scales::comma(sum(counts[str_detect(names(counts), "_1")]))
# general cohort
scales::comma(sum(counts[str_detect(names(counts), "_2")]))

table1_age_TRUE_REDACTED <- readr::read_csv(
  "release20220622/table1_age_TRUE_REDACTED.csv"
)

# percent of those aged 18-69 years in JCVI groups 4b & 6 in each cohort
# (5=cancer, 7=non-cancer)
table1_age_TRUE_REDACTED %>% 
  filter(Variable=="JCVI group", Characteristic %in% c("04b", "06")) %>% 
  arrange(Characteristic) %>%
  select(ends_with(as.character(c(5,7)))) %>%
  mutate(across(ends_with(as.character(c(5,7))),
                ~as.integer(str_remove(str_extract(.x, "\\d+\\%"), "\\%"))
                )
         ) %>%
  summarise(across(everything(), sum))

################################################################################
event_counts <- readr::read_csv(
  here::here(release_folder, "event_counts.csv")
)

event_counts %>%
  filter(
    outcome == "postest", 
    k %in% c(1,6),
    subgroup %in% c("cancer", "noncancer")
    ) %>%
  mutate(across(arm, ~if_else(arm == "unvax", .x, "vax"))) %>%
  group_by(subgroup, arm, k) %>%
  summarise(across(c("personyears", "n", "events"), ~sum(.x)), .groups = "keep") %>%
  ungroup() %>%
  select(-personyears) %>%
  rename(
    cohort = subgroup,
    group = arm, 
    period = k
  ) %>%
  mutate(
    absolute_risk = round(events/n, 5)#,
    # incidence_rate_per1000py = 1000*round(events/personyears,3)
    ) %>%
  kableExtra::kable("pipe")

event_counts %>%
  filter(
    outcome == "covidadmitted", 
    k %in% c(1,6),
    subgroup %in% c("cancer", "noncancer")
  ) %>%
  mutate(across(arm, ~if_else(arm == "unvax", .x, "vax"))) %>%
  group_by(subgroup, k) %>%
  summarise(across(c("personyears", "n", "events"), ~sum(.x)), .groups = "keep") %>%
  ungroup() %>%
  select(-personyears) %>%
  rename(
    cohort = subgroup,
    period = k
  ) %>%
  mutate(
    absolute_risk = round(events/n, 5)#,
    # incidence_rate_per1000py = 1000*round(events/personyears,3)
  ) %>%
  kableExtra::kable("pipe")

event_counts %>%
  filter(
    outcome == "covidadmitted", 
    k %in% c(1,6),
    subgroup %in% c("cancer", "noncancer") |
    subgroup %in% c("cancer_18-69", "cancer_70+", "noncancer_18-69", "noncancer_70+")
  ) %>%
  mutate(across(arm, ~if_else(arm == "unvax", .x, "vax"))) %>%
  group_by(subgroup, k) %>%
  summarise(across(c("personyears", "n", "events"), ~sum(.x)), .groups = "keep") %>%
  ungroup() %>%
  select(-personyears) %>%
  rename(
    cohort = subgroup,
    period = k
  ) %>%
  mutate(
    absolute_risk = round(events/n, 5)#,
    # incidence_rate_per1000py = 1000*round(events/personyears,3)
  ) %>%
  kableExtra::kable("pipe")



################################################################################
estimates_all <- readr::read_csv(
  here::here(release_folder, "estimates_all.csv")
) %>%
  mutate(across(subgroup, factor, 
                levels = seq_along(subgroups),
                labels = subgroups)) %>%
  filter(
    variable == "k",
    !reference_row,
    model=="max_adjusted",
    outcome %in% c("postest", "covidadmitted", "coviddeath", "noncoviddeath"),
    comparison=="both",
    prior
  ) %>%
  mutate(across(c(estimate, conf.low, conf.high), exp)) 

metareg_res <- readr::read_rds(here::here(release_folder, "metareg_res.rds")) %>%
  right_join(
    estimates_all %>%
      distinct(subgroup, comparison, outcome, model, prior)
  ) %>%
  select(subgroup, comparison, outcome, model, prior, starts_with("logrhr")) %>%
  mutate(across(starts_with("logrhr"), ~round(exp(.x),2)))

estimates_all %>%
  filter(
   period %in% c(1,6)
    ) %>%
  mutate(across(c(estimate, conf.low, conf.high), ~round(100*(1-.x),1)))  %>%
  transmute(
    subgroup, outcome, period, 
    value = glue("{estimate}% [{conf.high}, {conf.low}]")
  ) %>%
  pivot_wider(
    names_from = subgroup,
    values_from = value
  )


metareg_res %>%
  transmute(
    subgroup, outcome, 
    value = glue("{logrhr} [{logrhr_lower}, {logrhr_higher}]")
  ) %>%
  pivot_wider(
    names_from = subgroup,
    values_from = value
  )


# percentage point differences in 3-6 o 23-26
ppoint <- estimates_all %>%
  filter(
    as.integer(subgroup) %in% 1:2
    ) %>%
  mutate(across(c(estimate, conf.low, conf.high), ~round(100*(1-.x),1)))  

ppoint %>%
  mutate(text = glue("{format(estimate,nsmall=1,trim=TRUE)}% [{format(conf.high,nsmall=1,trim=TRUE)}, {format(conf.low,nsmall=1,trim=TRUE)}]")) %>%
  select(subgroup, outcome, period, estimate, text) %>%
  pivot_wider(
    names_from = period,
    values_from = c(estimate, text)
  ) %>%
  mutate(
    diff_1_6 = estimate_1 - estimate_6,
    diff_2_6 = estimate_2 - estimate_6
  ) %>%
  transmute(
    outcome,
    text_1_6 = glue("{format(diff_1_6,nsmall=1)} percentage points in our {subgroup} cohort ({text_1} during 3-6 weeks versus {text_6} during 23-26 weeks)"),
    text_2_6 = glue("{format(diff_2_6,nsmall=1)} percentage points in our {subgroup} cohort {text_2} during 3-6 weeks versus {text_6} during 23-26 weeks")
    ) %>%
  select(outcome, contains("1"))


###############################################################################
metareg_results_rhr <- readr::read_rds(
  here::here(release_folder, "metareg_results_rhr.rds")
)

metareg_results_rhr %>% 
  filter(
    subgroup%in% subgroups[1:2],
    comparison=="both", 
    outcome%in%c("postest", "covidadmitted", "coviddeath"),
    prior,
    model==1
  ) %>%
  mutate(across(starts_with("rhr"), ~round(.x,2))) %>%
  transmute(
    outcome, subgroup,
    value = glue("{rhr} [{rhr_lower}, {rhr_higher}]")
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
    subgroup%in% subgroups[5:8],
    comparison!="both", 
    outcome%in%c("postest","covidadmitted", "coviddeath"),
    prior,
    model==1
  ) %>%
  mutate(across(starts_with("rhr"), ~round(.x,2))) %>%
  mutate(value = glue("{rhr} [{rhr_lower}-{rhr_higher}]")) %>%
  select(comparison, outcome, subgroup, value) %>%
  pivot_wider(
    names_from = subgroup,
    values_from = value
  )
