################################################################################

# This script:
# applies the cox model for a given comparison and outcome


################################################################################
library(tidyverse)
library(glue)
library(survival)
library(broom)
library(broom.helpers)

## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  comparison <- "BNT162b2"
  subgroup_index <- 1L
  include_prior_infection <- TRUE
  outcome <- "anytest"
  
} else{
  comparison <- as.character(args[[1]])
  subgroup_index <- as.integer(args[[2]])
  include_prior_infection <- as.logical(args[[3]])
  outcome <- as.character(args[[4]])
}

################################################################################
# create directories
fs::dir_create(here::here("output", "models_cox", "data"))
fs::dir_create(here::here("output", "models_cox", "temp"))

################################################################################
# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))
subgroup <- subgroups[as.integer(str_extract(subgroup_index, "^\\d"))]

# read study parameters
study_parameters <- readr::read_rds(
  here::here("analysis", "lib", "study_parameters.rds"))
K <- study_parameters$K

################################################################################
model_glance <- list()
model_tidy <- list()

for (kk in 1:K) {
  
  # specify filename_suffix for saving models
  filename_suffix <- glue("{comparison}_{subgroup_index}_{include_prior_infection}_{outcome}_{kk}")
  
  # read data
  model_input <-  readr::read_rds(
    here::here("output", "preflight", "data", glue("model_input_{filename_suffix}.rds"))
  )
  cat(glue("k={kk}: is.null(model_input):"), "\n")
  print(is.null(model_input))
  
  
  
  
  if (is.null(model_input)) {
    # model input null if not enough events, see preflight.R for threshold
    # save empty tibbles to avoid errors
    
    model_glance[[kk]] <- tibble()
    model_tidy[[kk]] <- tibble()
    
    readr::write_rds(
      tibble(),
      here::here("output", "models_cox", "data", glue("modelcox_glance_{filename_suffix}.rds"))) 
    
    readr::write_rds(
      tibble(),
      here::here("output", "models_cox", "data", glue("modelcox_tidy_{filename_suffix}.rds"))) 
    
    readr::write_rds(
      tibble(),
      here::here("output", "models_cox", "data", glue("modelx_{comparison}_{filename_suffix}.rds")), 
      compress="gz")
    
    readr::write_rds(
      tibble(),
      here::here("output", "models_cox", "data", glue("model_object_empty_{filename_suffix}.rds")),
      compress = "gz"
    )
    
  } else {
    
    model_varlist <- readr::read_rds(
      here::here("analysis", "lib", "model_varlist.rds")
    )
    
    data_cox <- model_input$data
    formulas <- model_input$formulas
    
    ################################################################################
    # define formulas
    
    model_names = c(
      "1" = "Unadjusted",
      "2" = "Adjusted for age, sex and IMD",
      "3" = "Maximally adjusted"
    )
    
    formula_cox_1 <- formulas$unadjusted
    formula_cox_2 <- formulas$unadjusted %>%
      update(formulas$demographic)
    formula_cox_3 <- formula_cox_1 %>% 
      update(formulas$clinical)
    
    ################################################################################
    # apply and save models (coxph within functions was causing dramas with broom)
    
    # set max iterations
    opt_control <- coxph.control(iter.max = 30)
    
    coxmods <- list()
    cat(glue("...... fitting unadjusted model ......"), "\n")
    coxmods[[1]] <- try(coxph(
      formula = formula_cox_1,
      data = data_cox,
      # robust = TRUE,
      # id = patient_id,
      cluster = patient_id,
      na.action = "na.fail",
      control = opt_control))
    
    readr::write_rds(
      coxmods[[1]],
      here::here("output", "models_cox", "data", glue("model_object_1_{filename_suffix}.rds")),
      compress = "gz"
    )
    
    # cat(glue("...... adjusting for demographic vars ......"), "\n")
    coxmods[[2]] <- NULL
    #   try(coxph(
    #   formula = formula_cox_2,
    #   data = data_cox,
    #   # robust = TRUE,
    #   # id = patient_id,
    #   cluster = patient_id,
    #   na.action = "na.fail",
    #   control = opt_control))
    # 
    # readr::write_rds(
    #   coxmods[[2]],
    #   here::here("output", "models_cox", "data", glue("model_object_2_{filename_suffix}.rds")),
    #   compress = "gz"
    # )
    
    cat(glue("...... adjusting for demographic and clinical vars ......"), "\n")
    coxmods[[3]] <- try(coxph(
      formula = formula_cox_3,
      data = data_cox,
      # robust = TRUE,
      # id = patient_id,
      cluster = patient_id,
      na.action = "na.fail",
      control = opt_control))
    
    readr::write_rds(
      coxmods[[3]],
      here::here("output", "models_cox", "data", glue("model_object_3_{filename_suffix}.rds")),
      compress = "gz"
    )
    
    ################################################################################
    # process output
    glance <- list()
    tidy <- list()
    for (i in seq_along(coxmods)) {
      
      coxmod <- coxmods[[i]]
      
      if (!inherits(coxmod, "try-error") & !is.null(coxmod)) {
        
        glance[[i]] <-
          broom::glance(coxmod) %>%
          add_column(
            model = i,
            convergence = coxmod$info[["convergence"]],
            ram = format(object.size(coxmod), units="GB", standard="SI", digits=3L),
            .before = 1
          )
        
        tidy[[i]] <-
          broom.helpers::tidy_plus_plus(
            coxmod,
            exponentiate = FALSE
          ) %>%
          tidy_add_n(
            coxmod
          ) %>%
          add_column(
            model = i,
            .before=1
          ) %>%
          mutate(across(c(term, label),
                        ~if_else(
                          variable=="k" & term != "k0",
                          str_replace(.x, "\\d", as.character(kk)),
                          .x
                        )))
        
      } else {
        
        glance[[i]] <- tibble()
        tidy[[i]] <- tibble()
        
      }
      
    }
    
    model_glance_out <- bind_rows(glance) 
    readr::write_rds(
      model_glance_out,
      here::here("output", "models_cox", "data", glue("modelcox_glance_{filename_suffix}.rds"))) 
    
    model_glance[[kk]] <- model_glance_out %>%
      mutate(k = kk) %>%
      select(k, model, convergence) 
    
    model_tidy_out <- bind_rows(tidy)
    readr::write_rds(
      model_tidy_out,
      here::here("output", "models_cox", "data", glue("modelcox_tidy_{filename_suffix}.rds"))) 
    
    model_tidy[[kk]] <- model_tidy_out %>%
      filter(
        str_detect(term, "k\\d"),
        term != "k0"
      ) %>%
      select(model, term, estimate, std.error)
    
  }
  
}

# for quickly checking that models have run successfully for all k
model_glance <- bind_rows(model_glance)
if (!is_empty(model_glance)) {
  model_glance <- tibble(k = 1:K) %>% 
    left_join(model_glance, by = "k")
}
capture.output(
  model_glance %>% kableExtra::kable(format = "pipe"),
  file = here::here("output", "models_cox", "temp", glue("modelcox_glance_{comparison}_{subgroup_index}_{include_prior_infection}_{outcome}.txt")),
  append = FALSE
)

model_tidy <- bind_rows(model_tidy)
if (!is_empty(model_tidy)) {
  model_tidy <- tibble(term = str_c("k", 1:K)) %>% 
    left_join(model_tidy, by = "term")
}
capture.output(
  model_tidy %>% kableExtra::kable(format = "pipe"),
  file = here::here("output", "models_cox", "temp", glue("modelcox_tidy_{comparison}_{subgroup_index}_{include_prior_infection}_{outcome}.txt")),
  append = FALSE
)
