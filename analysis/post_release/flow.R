library(tidyverse)
library(glue)

################################################################################
# path for released results
release_folder <- here::here("release20220622")

################################################################################
# load data
eligibility_count_numeric <- readr::read_csv(
  here::here(release_folder, "eligibility_count_all.csv"))

################################################################################
# prepare data
n_groups <- eligibility_count_numeric %>%
  mutate(across(stage, ~str_replace(.x, "-", "_"))) %>%
  group_by(stage) %>%
  summarise(
    min = min(n, na.rm=TRUE), 
    max = max(n, na.rm=TRUE), 
    .groups = "keep"
    ) %>%
  ungroup() %>%
  pivot_wider(
    names_from = stage,
    values_from = c(min, max)
  )

eligibility_count <- eligibility_count_numeric %>%
  mutate(across(starts_with("n"),
                ~ str_c("(n = ", scales::comma(.x, accuracy=1), ")")))

################################################################################
# create tibble for results
ncol <- 5
nrow <- 100
flow <- matrix("", nrow=nrow, ncol=ncol)

# prior_covid <- str_c("(n = ", scales::comma(sum(eligibility_count_numeric[5:7,]$n_removed), accuracy = 1), ")")

# (meeting A but not meeting inclusion B or D) = 
# (meeting inclusion and exclusion A)
# -(meeting inclusion B)
# - (meeting inclusion D)
not_b_or_d <- n_groups$min_a_ex - n_groups$min_b_in - n_groups$min_d_in 
not_b_or_d <- str_c("(n =", scales::comma(not_b_or_d, accuracy = 1), ")")

# meeting A and B but not C
not_c <- n_groups$min_b_ex - n_groups$min_c_in
not_c <- str_c("(n =", scales::comma(not_c, accuracy = 1), ")")

vax_died_p1 <- scales::comma(eligibility_count_numeric[20,]$n_removed + eligibility_count_numeric[28,]$n_removed, accuracy = 1)
vax_dereg_p1 <- scales::comma(eligibility_count_numeric[21,]$n_removed + eligibility_count_numeric[29,]$n_removed, accuracy = 1)
vax_dose_p1 <- scales::comma(eligibility_count_numeric[22,]$n_removed + eligibility_count_numeric[30,]$n_removed, accuracy = 1)
unvax_died_p1 <- scales::comma(eligibility_count_numeric[14,]$n_removed + eligibility_count_numeric[24,]$n_removed, accuracy = 1)
unvax_dereg_p1 <- scales::comma(eligibility_count_numeric[15,]$n_removed + eligibility_count_numeric[25,]$n_removed, accuracy = 1)
unvax_dose_p1 <- scales::comma(eligibility_count_numeric[16,]$n_removed + eligibility_count_numeric[26,]$n_removed, accuracy = 1)

###
i <- 1
# A inc.
flow[i,3] <- glue("Meeting inclusion criteria A: Aged 18 years or over and registered for at least 1 year {eligibility_count[2,]$n}")
i <- i+2
# A ex.
flow[i,5] <- glue("Meeting exclusion criteria A: aged >120 {eligibility_count[3,]$n_removed}; missing  sex, ethnicity, region, or IMD {eligibility_count[4,]$n_removed}; in care home {eligibility_count[5,]$n_removed}; medically housebound {eligibility_count[6,]$n_removed}")
i <- i+2
# A remaining
flow[i,3] <- glue("Remaining: {eligibility_count[6,]$n}")
i <- i+2
###
# Not B or D
flow[i,5] <- glue("Not meeting inclusion criteria B or D: {not_b_or_d}")
i <- i+2
# B inc.
flow[i,2] <- glue("Meeting inclusion criteria B: second dose of BNT162b2 or ChAdOx received between 6 and 14 weeks after eligibility date {eligibility_count[7,]$n}")
i <- i+2
# B ex
flow[i,1] <- glue("Meeting exclusion criteria B: first dose received before eligibility date {eligibility_count[8,]$n_removed}; <6 or >=14 weeks between first and second dose {eligibility_count[9,]$n_removed}; healthcare worker {eligibility_count[10,]$n_removed}")
i <- i+2
###
# Not C
flow[i,1] <- glue("Not meeting inclusion criteria C: {not_c}") 
i <- i+2
# C inc.
flow[i,2] <- glue("Meeting inclusion criteria C: second dose received during SVP {eligibility_count[17,]$n}")
###
# D inc
flow[i,4] <- glue("Meeting inclusion criteria D: unvaccinated at start of SVP {eligibility_count[11,]$n}")
i <- i+2
###
# E ex vax
flow[i,1] <- glue("Meeting exclusion citeria E: resident in care home before start of SVP {eligibility_count[18,]$n_removed}; end of life care before start of SVP {eligibility_count[19,]$n_removed}")
# E ex unvax
flow[i,5] <- glue("Meeting exclusion citeria E: resident in care home before start of SVP {eligibility_count[12,]$n_removed}; end of life care before start of SVP {eligibility_count[13,]$n_removed}")
i <- i+2
###
# P1 exclusions vax
flow[i,1] <- glue("Occurrence of the following before the start of period 1: death {vax_died_p1}; de-registration {vax_dereg_p1}; subsequent dose {vax_dose_p1}")
# P1 exclusions unvax
flow[i,5] <- glue("Occurrence of the following before the start of period 1: death {unvax_died_p1}; de-registration {unvax_dereg_p1}; subsequent dose {unvax_dose_p1}")
# remaining
i <- i+2
flow[i,2] <- glue("Eligible for vaccinated group in period 1 {eligibility_count[30,]$n}")
flow[i,4] <- glue("Eligible for unvaccinated group in period 1 {eligibility_count[26,]$n}")

flow <- flow[1:i,1:5]

readr::write_delim(
  as.data.frame(flow),
  file = file.path(release_folder, "flow.csv"),
  delim = "\t"
)

