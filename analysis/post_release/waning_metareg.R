# metaregression

library(tidyverse)

release_folder <- here::here("release20220622")

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
data_metareg <- estimates_all %>%
  filter(
    outcome %in% outcomes,
    !reference_row,
    variable %in% "k"
  ) %>%
  mutate(
    subgroup = as.integer(str_extract(subgroup, "\\d")),
    seloghr = (conf.high - conf.low)/(2*qnorm(0.975)),
    k = as.integer(label) - 1 # so that intercept at comparison period 1
  ) %>%
  select(subgroup, comparison, outcome, model, prior, k, loghr = estimate, seloghr) 

################################################################################

data_metareg_split <- data_metareg %>%
  group_split(subgroup, comparison, outcome, model, prior)


my_rma <- function(.data) {
  
  data_out <- .data %>%
    distinct(subgroup, comparison, outcome, model, prior) 
    
  # check >2 estimates
  if(length(.data$seloghr) > 2) {
    metafor_res <- try(metafor::rma(
      yi = loghr,
      sei = seloghr,
      mods = ~k,
      method = "REML", # restricted maximum likelihood estimator (Viechtbauer, 2005; Raudenbush, 2009)
      test = "knha", # Knapp-Hartung method
      data = .data
    ))
    
    data_out <- data_out %>%
      bind_cols(
        tibble(
          logrhr = metafor_res$beta[2,1],
          selogrhr = metafor_res$se[2],
          loghr1 = metafor_res$beta[1,1],
          seloghr1 = metafor_res$se[1],
          tau2 = metafor_res$tau2,
          I2 = metafor_res$I2,
          H2 = metafor_res$H2,
          R2 = metafor_res$R2
        )
      )
  } 
  
  return(data_out)
  
}

rma_res <- bind_rows(lapply(
  data_metareg_split,
  my_rma
)) 

metareg_res <- rma_res %>%
  mutate(across(subgroup,
                factor,
                levels = seq_along(subgroups),
                labels = subgroups)) %>%
  mutate(across(comparison, 
                factor,
                levels = c("BNT162b2", "ChAdOx1", "both")
  )) %>%
  mutate(across(outcome,
                factor,
                levels = c("covidadmitted",
                           "coviddeath",
                           "postest",
                           "noncoviddeath",
                           "anytest"))) %>%
  mutate(
    # slope ci
    logrhr_lower = logrhr - qnorm(0.975)*selogrhr,
    logrhr_higher = logrhr + qnorm(0.975)*selogrhr,
    # intercept ci
    loghr1_lower = loghr1 - qnorm(0.975)*seloghr1,
    loghr1_higher = loghr1 + qnorm(0.975)*seloghr1,
  ) 

readr::write_rds(
  metareg_res,
  here::here(release_folder, "metareg_res.rds")
)

