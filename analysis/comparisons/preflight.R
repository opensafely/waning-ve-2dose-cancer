################################################################################

# This script:


################################################################################
library(tidyverse)
library(glue)
library(gt)

## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  comparison <- "both"
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
# create output directories
fs::dir_create(here::here("output", "preflight", "data"))
fs::dir_create(here::here("output", "preflight", "tables"))

################################################################################
# read study parameters
study_parameters <- readr::read_rds(
  here::here("analysis", "lib", "study_parameters.rds"))
K <- study_parameters$K

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))

# model covariates
model_varlist <- readr::read_rds(
  here::here("analysis", "lib", "model_varlist.rds")
)
vars <- unname(c(model_varlist$demographic, model_varlist$clinical))

################################################################################
## read data
# covariates data
data_covs <- readr::read_rds(
  here::here("output", "data", "data_all.rds")) %>%
  select(patient_id, 
         jcvi_group, elig_date, region,
         unname(model_varlist$demographic),
         unname(model_varlist$clinical))

# time to event data
data_tte_0 <- readr::read_rds(
  here::here("output", "tte", "data", glue("data_tte_{outcome}.rds"))) 

# select subgroup
if (subgroup_index == 1) {
  data_tte_1 <- data_tte_0 %>%
    filter(
      cancer_subgroup != "noncancer"
    ) 
} else if (subgroup_index == 2) {
  data_tte_1 <- data_tte_0 %>%
    filter(
      cancer_subgroup == "noncancer"
    ) 
} else if (subgroup_index == 3) {
  data_tte_1 <- data_tte_0 %>%
    filter(
      cancer_subgroup == "cancer_haem"
    ) 
} else if (subgroup_index == 4) {
  data_tte_1 <- data_tte_0 %>%
    filter(
      cancer_subgroup == "cancer_solid"
    ) 
} else if (subgroup_index == 5) {
  data_tte_1 <- data_tte_0 %>%
    filter(
      cancer_subgroup != "noncancer",
      age_subgroup == "18-69"
    ) 
} else if (subgroup_index == 6) {
  data_tte_1 <- data_tte_0 %>%
    filter(
      cancer_subgroup != "noncancer",
      age_subgroup == "70+"
    ) 
} else if (subgroup_index == 7) {
  data_tte_1 <- data_tte_0 %>%
    filter(
      cancer_subgroup == "noncancer",
      age_subgroup == "18-69"
    ) 
} else if (subgroup_index == 8) {
  data_tte_1 <- data_tte_0 %>%
    filter(
      cancer_subgroup == "noncancer",
      age_subgroup == "70+"
    ) 
}

# select the arms to include
# derive weights so that BNT162b2 and ChAdOx1 equal weights for "both" analysis
if (comparison == "both") {
  
  # derive weights for each brand
  # weights for smaller group = 1
  # weights for larger group = n_smaller group / n_larger group
  arm_weights <- data_tte_1 %>%
    filter(arm != "unvax") %>%
    group_by(arm) %>%
    count() %>%
    ungroup() %>%
    mutate(weight = min(n)/n)
  
  # grouping brands together
  data_tte_2 <- data_tte_1 %>%
    left_join(arm_weights, by = "arm") %>%
    mutate(across(weight, ~if_else(arm == "unvax", 1, .x))) %>%
    mutate(across(arm, 
                  ~factor(
                    if_else(.x == "unvax", as.character(.x), "vax"),
                    levels = c("vax", "unvax")
                  ))) 
    
    
} else {
  # single brand
  data_tte_2 <- data_tte_1 %>%
    filter(arm %in% c("unvax", comparison)) %>%
    mutate(weight = 1)
}

# apply prior infection criteria
if (include_prior_infection) {
  data_tte_3 <- data_tte_2
} else {
  data_tte_3 <- data_tte_2 %>%
    filter(!prior_infection_subgroup)
}

################################################################################
# read functions

# redaction functions
source(here::here("analysis", "functions", "redaction_functions.R"))

################################################################################
# join tte data to covariates
data_0 <- data_tte_3 %>%
  select(-ends_with("_subgroup")) %>%
  left_join(data_covs, 
            by = "patient_id") %>%
  mutate(strata_var = factor(str_c(jcvi_group, elig_date, region, sep = ", "))) %>%
  droplevels()

events_threshold <- 2
# reasons for applying the events threshold are twofold:
# 1. avoids issues with statistical convergence
# 2. ensuring that >2 events in each cell of the cross-tab of 
#    exposure variable x categorical covariate also ensures that >5 events per
#    level of each covariate


# keep only periods with > 2 events in each arm
events_per_period <- data_0 %>%
  group_by(k, arm) %>%
  summarise(events = sum(status), .groups="keep") %>%
  ungroup(arm) %>%
  summarise(
    n = n(), # number of arms with events
    min_events = min(events), # min event across arms
    .groups = "keep") %>%
  ungroup() %>%
  mutate(keep = min_events > events_threshold)

keep_periods <- as.integer(events_per_period$k[events_per_period$keep])
drop_periods <- as.integer(events_per_period$k[!events_per_period$keep])

data_00 <- data_0 %>%
  filter(k %in% keep_periods) %>%
  droplevels() 

# prepare data for each period k
for (kk in 1:K) {
  
  data_1 <- data_00 %>%
    filter(k == kk)
  
  # do not run for kk if all periods dropped
  if (nrow(data_1) == 0) {
    readr::write_file(
      x="",
      here::here("output", "preflight", "tables", glue("eventcheck_{comparison}_{subgroup_index}_{include_prior_infection}_{outcome}_{kk}.html")),
      append = FALSE
    )
    readr::write_rds(
      NULL,
      here::here("output", "preflight", "data", glue("model_input_{comparison}_{subgroup_index}_{include_prior_infection}_{outcome}_{kk}.rds"))
    )
    capture.output(
      print("No events"),
      file = here::here("output", "preflight", "tables", glue("preflight_report_{comparison}_{subgroup_index}_{include_prior_infection}_{outcome}_{kk}.txt")),
      append = FALSE
    )
    next
  }
    
    # only keep categorical covariates with > 2 events per level per arm
    
    ################################################################################
    # tabulate events per level
    cat("...split data by k and status...\n")
    tbl_list <- data_1 %>%
      select(k, status, all_of(vars)) %>%
      select(-age) %>%
      group_split(status)
    
    # names for each element in list
    group_split_labels <- lapply(
      tbl_list,
      function(x) str_c(unique(x$k), unique(x$status))
    ) %>% 
      unlist()
    
    cat("...summarise number of events...\n")
    # summarise the number of events by level of covariates (within comparisons)
    tbltab_list <- tbl_list %>%
      map(~.[,-c(1,2)]) %>%
      map(
        function(data){
          map(data,
              function(x) {
                tab <- table(x)
                tibble(.level = names(tab),
                       n = as.vector(tab))
              }) %>%
            bind_rows(.id="variable")
        }
      )
    
    # apply names
    names(tbltab_list) <- group_split_labels
    
    cat("...prepare table...\n")
    tbltab <- bind_rows(
      tbltab_list,
      .id = "group"
    ) %>%
      mutate(
        k = str_extract(group, "\\d"),
        status = as.logical(str_remove(group, "\\d")))  %>%
      select(-group) %>%
      pivot_wider(
        names_from = status,
        values_from = n
      ) %>%
      pivot_wider(
        names_from = k,
        values_from = c("FALSE", "TRUE"),
        names_glue = "period{k}_{.value}"
      ) %>%
      mutate(across(starts_with("period"), 
                    ~ if_else(is.na(.x), 0L, .x))) %>% 
      group_by(variable) %>%
      mutate(across(starts_with("period"), redactor2)) %>% 
      mutate(across(starts_with("period"), 
                    ~ if_else(is.na(.x), "-", scales::comma(.x, accuracy = 1)))) %>% 
      ungroup()
    
    cat("...format and save table...\n")
    tbltab %>%
      gt(
        groupname_col="variable",
        rowname_col = ".level"
      ) %>%
      tab_spanner_delim("_") %>%
      tab_stubhead(label = "variable") %>%
      opt_css(css = ".gt_stub { padding-left: 50px !important; }") %>%
      tab_style(
        style = list(
          cell_fill(color = "lightcyan")
        ),
        locations = cells_body(
          columns = ends_with("TRUE")
        )
      ) %>%
      tab_style(
        style = list(
          cell_fill(color = "lightcyan")
        ),
        locations = cells_column_labels(
          columns = ends_with("TRUE")
        )
      ) %>%
      gtsave(
        filename = glue("eventcheck_{comparison}_{subgroup_index}_{include_prior_infection}_{outcome}_{kk}_REDACTED.html"),
        path = here::here("output", "preflight", "tables")
      )
    
    
    ################################################################################
    # check levels per covariate
    n_levels <- sapply(
      data_1 %>% 
        select(all_of(vars)) %>% 
        select(-age) %>%
        mutate(across(everything(), as.factor)),
      function(x) length(levels(x))
    )
    
    # calculate events per level during each period, split by treatment arm
    data_1_list <- data_1 %>% 
      arrange(k, arm) %>%
      group_split(k, arm) %>%
      as.list()
    
    events_per_level_list <- list()
    for (i in seq_along(data_1_list)) {
      
      events_per_level_list[[i]] <- data_1_list[[i]] %>%
        filter(status) %>%
        select(all_of(vars)) %>%
        select(-age) %>%
        map(function(x) {
          tab <- table(x)
          tibble(level = names(tab), 
                 n = as.vector(tab))
        }) %>%
        bind_rows(.id="variable") %>%
        mutate(keep = n > events_threshold) %>%
        group_by(variable) %>%
        mutate(
          n_levels = n(),
          n_keep = sum(keep)
        ) %>%
        ungroup() 
    }
    
    drop_vars_list <- list()
    for (i in seq_along(data_1_list)) {
      
      drop_vars_list[[i]] <- events_per_level_list[[i]] %>%
        filter(
          # remove if one level or binary and one level with too few events
          # would be one level only if other levels dropped previously with droplevels
          n_levels ==1 | (n_levels == 2  & n_keep == 1) 
        ) %>% 
        distinct(variable)
      
    }
    drop_vars <- bind_rows(drop_vars_list) %>% distinct() %>% unlist() %>% unname()
    
    ################################################################################
    merge_levels <- function(var) {
      
      data <- data_1 %>%
        select(status, k, arm, all_of(var))
      
      event_count_fun <- function(var) {
        
        data %>%
          filter(status) %>%
          group_by(k, arm, !! sym(var)) %>%
          count() %>%
          ungroup() %>%
          summarise(min_n=min(n)) %>%
          unlist() %>%
          unname()
        
      }
      
      var_levs <- levels(data[[var]])
      
      for (i in rev(seq_along(var_levs))) {
        
        if (event_count_fun(var) > events_threshold) {
          # no merging required if minimum count is above threshold
          break
        } else if (i > 2) {
          # merge labels i and i-1
          var_levs[(i-1)] <- str_c(var_levs[(i-1)], var_levs[i], sep = " / ")
          var_levs <- var_levs[-i]
          # merge levels in data  
          data <- data %>%
            mutate(across(var,
                          ~ factor(
                            if_else(
                              as.integer(.x) == i,
                              as.integer(i-1),
                              as.integer(.x)),
                            levels = 1:(i-1),
                            labels = var_levs)))
          
          
        } else {
          # if loop reaches i=2, add var to drop_vars and break out of loop
          drop_vars <- c(drop_vars, var)
          break
        }
      }
      data %>% select(all_of(var))
    }
    
    merged_vars <- lapply(names(n_levels[n_levels>2]), merge_levels)
    
    data_2 <- data_1 %>%
      select(-all_of(names(n_levels[n_levels>2]))) %>%
      bind_cols(merged_vars) %>%
      select(-all_of(drop_vars)) %>%
      droplevels()
    
    ################################################################################
    # create exposure variable
    data_3 <- data_2 %>%
      mutate(across(k,
                    ~factor(
                      if_else(arm %in% "unvax",
                              0L,
                              as.integer(k)))))
    
    ################################################################################
    # create age variables
    data_5 <- data_3 %>% mutate(age2=age^2) 
    
    ################################################################################
    # define formulas
    
    formula_unadj <- formula(Surv(tstart, tstop, status, type = "counting") ~ k + strata(strata_var))
    
    demog_vars <- unname(model_varlist$demographic[which(model_varlist$demographic %in% names(data_5))])
    demog_vars <- demog_vars[demog_vars != "age"]
    demog_vars <- c("age", "age2", demog_vars)
    formula_demog <- as.formula(str_c(c(". ~ . ", demog_vars), collapse = " + "))
    
    clinical_vars <- unname(model_varlist$clinical[which(model_varlist$clinical %in% names(data_5))])
    formula_clinical <- as.formula(str_c(c(". ~ . ", clinical_vars), collapse = " + "))
    
    ################################################################################
    
    formulas_list <-  list(
      "unadjusted" = formula_unadj, 
      "demographic" = formula_demog, 
      "clinical" = formula_clinical)
    
    model_input <- list(
      data = data_5,
      formulas = formulas_list
    )
    
    readr::write_rds(
      model_input,
      here::here("output", "preflight", "data", glue("model_input_{comparison}_{subgroup_index}_{include_prior_infection}_{outcome}_{kk}.rds"))
    )
    
    ################################################################################
    
    preflight_report <- function(
      dropped_variables,
      merged_variables,
      subgroup_string = subgroups[subgroup_index]
    ) {
      ####
      cat(glue("Comparison = {comparison}; Subgroup = {subgroup_string}; Include prior covid = {include_prior_infection}; Outcome = {outcome}; Comparison period = {kk}"), "\n")
      ####
      cat("---\n")
      if (is_empty(dropped_variables)) {
        dropped_variables <- "none"
      } else {
        dropped_variables <- str_c(str_c("- ", dropped_variables), collapse = "\n")
      }
      cat(glue("Dropped variables:\n{dropped_variables}"), "\n")
      ####
      cat("---\n")
      cat("Ordinal variable categories:", "\n")
      levs <- bind_cols(merged_variables) %>%
        select(-all_of(dropped_variables[dropped_variables %in% names(bind_cols(merged_variables))])) %>%
        map(levels)
      levs_tidy <- lapply(names(levs),
                          function(x) {
                            cat("\n---", x, "---\n")
                            cat(str_c(levs[[x]], collapse = "\n")) 
                            cat("\n")
                          })
      ####
      cat("\n")
      cat("---\n")
      cat(glue("Formulas:"), "\n")
      print(formulas_list)
    }
    
    capture.output(
      preflight_report(
        dropped_variables = drop_vars,
        merged_variables = merged_vars
      ),
      file = here::here("output", "preflight", "tables", glue("preflight_report_{comparison}_{subgroup_index}_{include_prior_infection}_{outcome}_{kk}.txt")),
      append = FALSE
    )
    
  
} 
