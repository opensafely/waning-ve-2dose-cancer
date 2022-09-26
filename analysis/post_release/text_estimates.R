library(tidyverse)

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
# extract group x cohort counts
counts_formatted <- unlist(table1_main_TRUE_REDACTED[1,3:8])
# remove commas and reformat as numeric
counts <- as.numeric(str_remove_all(counts_formatted, ","))
names(counts) <- cols

## number eligibile for each vaccine group
for (i in c("BNT162b2", "ChAdOx1", "Unvaccinated")) {
  total <- sum(counts[str_detect(names(counts), glue::glue("^{i}"))])
  pct_cancer <- round(100*counts[glue::glue("{i}_1")]/total,0)
  cat(
    glue::glue("{i}:"), "\n",
    glue::glue("N eligible = {scales::comma(total, accuracy=1)}"), "\n",
    glue::glue("% cancer = {pct_cancer}"), "\n"
    )
}

## total in each cohort
# cancer cohort:
scales::comma(sum(counts[str_detect(names(counts), "_1")]))
# general cohort
scales::comma(sum(counts[str_detect(names(counts), "_2")]))

# read summary stats in the cancer and general cohorts split by age subgroup
table1_age_TRUE_REDACTED <- readr::read_csv(
  "release20220622/table1_age_TRUE_REDACTED.csv"
)

# percent of those aged 18-69 years in JCVI groups 4b & 6 in each cohort
# (5=cancer, 7=non-cancer)
table1_age_TRUE_REDACTED %>% 
  # keep rows for groups 4b & 6
  filter(Variable=="JCVI group", Characteristic %in% c("04b", "06")) %>% 
  arrange(Characteristic) %>%
  # keep 18-69 subgroups
  select(ends_with(as.character(c(5,7)))) %>%
  # remove formatting and code as integers
  mutate(across(ends_with(as.character(c(5,7))),
                ~as.integer(str_remove(str_extract(.x, "\\d+\\%"), "\\%"))
                )
         ) %>%
  # sum percentages across the JCVI groups
  summarise(across(everything(), sum))

################################################################################
# read event counts
event_counts <- readr::read_csv(
  here::here(release_folder, "event_counts.csv")
)

# calculate absolute risk in comparison periods 1 and 6
# event_counts %>%
#   filter(
#     outcome == "covidadmitted", 
#     k %in% c(1,6),
#     subgroup %in% c("cancer", "noncancer") |
#     subgroup %in% c("cancer_18-69", "cancer_70+", "noncancer_18-69", "noncancer_70+")
#   ) %>%
#   mutate(across(arm, ~if_else(arm == "unvax", .x, "vax"))) %>%
#   group_by(subgroup, k) %>%
#   summarise(across(c("personyears", "n", "events"), ~sum(.x)), .groups = "keep") %>%
#   ungroup() %>%
#   select(-personyears) %>%
#   rename(
#     cohort = subgroup,
#     period = k
#   ) %>%
#   mutate(
#     absolute_risk = round(events/n, 5)#,
#     # incidence_rate_per1000py = 1000*round(events/personyears,3)
#   ) %>%
#   kableExtra::kable("pipe")

################################################################################
# read all hr estimates and filter estimates of vaccine effectiveness
estimates_all <- readr::read_csv(
  here::here(release_folder, "estimates_all.csv")
) %>%
  mutate(across(subgroup, factor, 
                levels = seq_along(subgroups),
                labels = subgroups)) %>%
  filter(
    variable == "k", # comparison period estimates
    !reference_row, 
    model=="max_adjusted", # keep only max adjusted model
    outcome %in% c("postest", "covidadmitted", "coviddeath", "noncoviddeath"),
    comparison=="both", # keep only combined effectiveness for BNT162b2 and ChAdOx1
    prior # include those with prior infection
  ) %>%
  # calculate exponential of estimate and confidence intervals
  mutate(across(c(estimate, conf.low, conf.high), exp)) 


# percentage point differences from 3-6 to 23-26 
ppoint <- estimates_all %>%
  filter(
    as.integer(subgroup) %in% 1:2 # cancer and general cohorts
    ) %>%
  mutate(across(c(estimate, conf.low, conf.high), ~round(100*(1-.x),1)))  

# generate text to paste into manuscript
ppoint %>%
  mutate(text = glue::glue("{format(estimate,nsmall=1,trim=TRUE)}% [{format(conf.high,nsmall=1,trim=TRUE)}, {format(conf.low,nsmall=1,trim=TRUE)}]")) %>%
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
    text_1_6 = glue::glue("{format(diff_1_6,nsmall=1)} percentage points in our {subgroup} cohort ({text_1} during 3-6 weeks versus {text_6} during 23-26 weeks)"),
    text_2_6 = glue::glue("{format(diff_2_6,nsmall=1)} percentage points in our {subgroup} cohort {text_2} during 3-6 weeks versus {text_6} during 23-26 weeks")
    ) %>%
  select(outcome, contains("1"))


###############################################################################
# read ratio of HRs data
metareg_results_rhr <- readr::read_rds(
  here::here(release_folder, "metareg_results_rhr.rds")
)

print_rhrs <- function(i) {
  metareg_results_rhr %>% 
    filter(
      subgroup%in% subgroups[i],
      comparison=="both", 
      outcome%in%c("postest", "covidadmitted", "coviddeath"),
      prior,
      model==3
    ) %>%
    mutate(across(starts_with("rhr"), ~round(.x,2))) %>%
    transmute(
      outcome, subgroup,
      value = glue::glue("{rhr} [{rhr_lower}, {rhr_higher}]")
    ) 
}

# ratios of HRs in cancer and general cohorts
print_rhrs(1:2)

# ratios of HRs in hame and solid subgroups
print_rhrs(3:4)

# ratios of HRs in cancer and general cohorts and age subgroups
print_rhrs(5:8)
