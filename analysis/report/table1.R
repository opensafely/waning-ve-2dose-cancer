library(tidyverse)
library(glue)
library(lubridate)
library(gt)

################################################################################
fs::dir_create(here::here("output", "tables"))
fs::dir_create(here::here("output", "report", "tables"))

################################################################################
# read list of covariates for model
model_varlist <- readr::read_rds(
  here::here("analysis", "lib", "model_varlist.rds"))

# read comparisons
comparisons <- readr::read_rds(
  here::here("analysis", "lib", "comparisons.rds"))

# read strata_vars
strata_vars <- readr::read_rds(
  here::here("analysis", "lib", "strata_vars.rds"))
strata_vars <- strata_vars[strata_vars!="elig_date"]

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))

################################################################################
# processed data
data_all <- readr::read_rds(
  here::here("output", "data", "data_all.rds")) 

################################################################################
# redaction functions
source(here::here("analysis", "functions", "redaction_functions.R"))

# function to be applied in dplyr::filter
no_evidence_of <- function(cov_date, index_date) {
  is.na(cov_date) | index_date < cov_date
}

################################################################################
# prepare data and create eligibility count tibble for flow diagram
data_tables_0 <- data_all %>%
  mutate(group = if_else(arm == "unvax", "unvax", "vax"))

# eligible for comparison period 1
eligibility_count <- data_tables_0 %>% 
  group_by(group) %>%
  count() %>%
  ungroup() %>%
  transmute(
    description = glue("{group}: satisfying eligibility criteria up to and including box E."),
    n
  )

# remove if death before start of comparison 1 (includes day 0)
data_tables <- data_tables_0 %>%
  filter_at(
    vars("death_date"),
    all_vars(no_evidence_of(., start_1_date))) 

eligibility_count <- eligibility_count %>%
  bind_rows(
    data_tables %>% 
      group_by(group) %>%
      count() %>%
      ungroup() %>%
      transmute(
        description = glue("{group}: after removing those who died before start of period 1."),
        n
      ))

# remove if dereg before start of comparison 1
data_tables <- data_tables %>%
  filter_at(
    vars("dereg_date"),
    all_vars(no_evidence_of(., start_1_date))) 

eligibility_count <- eligibility_count %>%
  bind_rows(
    data_tables %>% 
      group_by(group) %>%
      count() %>%
      ungroup() %>%
      transmute(
        description = glue("{group}: after removing those who deregistered before start of period 1."),
        n
      ))

# remove if subsequent_vax before start of comparison 1
data_tables <- data_tables %>%
  filter_at(
    vars("subsequent_vax_date"),
    all_vars(no_evidence_of(., start_1_date))) 

eligibility_count <- eligibility_count %>%
  bind_rows(
    data_tables %>% 
      group_by(group) %>%
      count() %>%
      ungroup() %>%
      transmute(
        description = glue("{group}: after removing those who received a subsequent dose before start of period 1."),
        n
      ))

eligibility_count_p1 <- eligibility_count %>%
  mutate(group = str_extract(description, "\\w+:")) %>%
  arrange(group) %>%
  # round up to nearest 7
  mutate(across(n, ~ceiling_any(.x, to=7))) %>%
  group_by(group) %>%
  mutate(n_removed = lag(n) - n) %>%
  ungroup() %>%
  select(-group)

readr::write_csv(
  eligibility_count_p1,
  here::here("output", "tables", "eligibility_count_p1.csv"))

################################################################################
# combine and save eligibility count tables
eligibility_count_all <- bind_rows(
  readr::read_csv(here::here("output", "tables", "eligibility_count_ab.csv")),
  readr::read_csv(here::here("output", "tables", "eligibility_count_cde.csv")),
  eligibility_count_p1 %>% mutate(stage="p1")
)

readr::write_csv(
  eligibility_count_all,
  here::here("output", "tables", "eligibility_count_all.csv"))

################################################################################
# function for table1 for each analysis
make_table1 <- function(data, include_prior_infection = TRUE) {
  
  # get the name of the data object to label the analysis
  data_name <- str_remove(deparse(substitute(data)), "data_")
  
  cat(glue("---- {data_name} ----\n"))
  
  # remove those with prior infection if necessary
  if (!include_prior_infection) {
    data <- data %>% filter(!prior_infection_subgroup)
  }
  
  cat("---- split dataset by subgroup ----\n")
  data <- data %>%
    select(patient_id, arm, region, jcvi_group, subgroup, ethnicity,
           all_of(c(unname(model_varlist$demographic), unname(model_varlist$clinical)))) %>%
    group_split(subgroup)
  
  
  cat("---- define obejcts ----\n")
  variables <- c(unname(strata_vars), unname(c(model_varlist$demographic, "ethnicity", model_varlist$clinical)))
  variables <- variables[variables != "age"]
  vars_ordered_levs <- c("region", "jcvi_group", "sex", "imd", "ethnicity", "bmi", "multimorb", "test_hist_n")
  
  
  ethnicity_var <- "ethnicity"
  names(ethnicity_var) <- "Ethnicity"
  # tibble for assigning tidy variable names
  var_labels <- tibble(
    variable = c(
      strata_vars,
      model_varlist$clinical, 
      # model_varlist$multimorb, 
      ethnicity_var,
      model_varlist$demographic),
    variable_label = names(c(
      strata_vars, 
      model_varlist$clinical, 
      # model_varlist$multimorb,
      ethnicity_var,
      model_varlist$demographic))
  )
  
  # define function to summarise each variable and (within variables) 
  summary_var <- function(.data, var) {
    out <- .data %>%
      group_by(arm, !! sym(var)) %>%
      count() %>%
      ungroup(!! sym(var)) %>%
      mutate(arm_total = sum(n)) %>%
      ungroup() %>%
      # round all frequencies up to nearest 7
      mutate(across(n, ~ceiling_any(.x, to=7))) %>%
      mutate(percent = round(100*n/arm_total,0)) %>%
      mutate(across(n, ~scales::comma(.x, accuracy = 1))) %>%
      mutate(value = as.character(glue("{n} ({percent}%)"))) %>%
      select(arm, !! sym(var), value) %>%
      pivot_wider(
        names_from = arm, 
        values_from = value
      ) %>%
      mutate(variable = var) %>%
      rename("category" = !! sym(var))
    
    if (is.logical(out$category)) {
      out <- out %>%
        filter(category) %>%
        mutate(across(category, ~ "yes"))
    }
    
    out %>% mutate(across(category, as.character))
    
  }
  
  # empty list for output
  table1_tidy_n <- list()
  # subgroups in analysis
  table1_subgroups <- sapply(data, function(x) x$subgroup[[1]])
  # for each subgroup in the data
  for (i in seq_along(data)) {
    # capture subgroup label
    subgroup_labels <- data[[i]]$subgroup[[1]]
   
    # define function for creating tibble of categories for each variable
    var_tibble <- function(var) {
      var_region <- tibble(
        variable = var,
        category = levels(data[[i]][[var]])
      )
    }
    
    # tibble for specifying order of variables and categories
    var_order <- tibble(
      variable = variables
    ) %>%
      left_join(
        bind_rows(
          lapply(vars_ordered_levs,
                 var_tibble)),
        by = "variable"
      ) %>%
      mutate(across(category, ~ if_else(is.na(.x), "yes", .x)))
    
    cat("---- make table 1 ----\n")
    # make table1
    table1 <- bind_rows(lapply(
      variables,
      function(x)
        data[[i]] %>% summary_var(var = x)
    ))
    
    cat("---- tidy table 1 ----\n")
    table1_tidy <- var_order %>% 
      left_join(var_labels, by = "variable") %>%
      left_join(table1, by = c("category", "variable")) %>%
      rename(Variable = variable_label, Characteristic = category, Unvaccinated = unvax) %>%
      select(-variable) %>%
      select(Variable, Characteristic, everything()) %>%
      mutate(across(c(BNT162b2, ChAdOx1, Unvaccinated), 
                    ~ if_else(is.na(.x), "- (-%)", .x))) 
    
    cat("---- age summary ----\n")
    age_summary <- data[[i]] %>%
      group_by(arm) %>%
      summarise(
        median = median(age, na.rm=TRUE),
        iqr = IQR(age, na.rm = TRUE),
        .groups = "keep") %>% 
      ungroup() %>%
      transmute(arm, value = as.character(glue("{median} ({iqr})"))) %>%
      pivot_wider(names_from = "arm", values_from = "value") %>%
      mutate(Variable = "Age", Characteristic = "Median (IQR)") 
    
    cat("---- bind table1 parts ----\n")
    table1_tidy_n[[i]] <- data[[i]] %>% 
      group_by(arm) %>% 
      count() %>% 
      ungroup() %>% 
      # round total counts up to nearest 7
      mutate(across(n, ~ceiling_any(.x, to=7))) %>%
      pivot_wider(names_from = arm, values_from = n) %>% 
      mutate(Variable = "", Characteristic = "N") %>%
      mutate(across(c(BNT162b2, ChAdOx1, unvax), 
                    ~ scales::comma(.x, accuracy = 1))) %>%
      bind_rows(
        age_summary
      ) %>%
      rename(Unvaccinated = unvax) %>%
      bind_rows(
        table1_tidy
      )
     
  }
  
  cat("---- join table1 subgroups ----\n")
  out <- table1_tidy_n[[1]] %>%
    rename_with(
      .fn = ~str_c(.x, "_", which(subgroups == table1_subgroups[1])),
      .cols = any_of(c(comparisons, "Unvaccinated"))
      )
  for (i in 2:length(data)) {
    out <- out %>% 
      left_join(
        table1_tidy_n[[i]] %>%
          rename_with(
            .fn = ~str_c(.x, "_", which(subgroups == table1_subgroups[i])),
            .cols = any_of(c(comparisons, "Unvaccinated"))
            ),
            .cols = any_of(c(comparisons, "Unvaccinated")), 
        by = c("Variable", "Characteristic")
        )
  }
  
  out <- out %>% 
    select(
      Variable, Characteristic,
      ends_with(as.character(which(subgroups %in% table1_subgroups)))
    )
  
  cat("---- save table1.csv ----\n")
  # save table1_tidy
  readr::write_csv(
    out,
    here::here("output", "report", "tables", glue("table1_{data_name}_{include_prior_infection}_REDACTED.csv"))
    )
  
  cat("---- save table1.html ----\n")
  out %>%
    gt(
      groupname_col="Variable",
      rowname_col = "Characteristic"
    ) %>%
    tab_header(
      title = glue("Analysis: {data_name}"),
      subtitle = "Patient characteristics as of second vaccination period + 2 weeks") %>%
    tab_style(
      style = cell_text(weight="bold"),
      locations = list(
        cells_column_labels(
          columns = everything()
        ),
        cells_row_groups(
          groups = everything()
        ))
    ) %>%
    gtsave(
      filename = glue("table1_{data_name}_{include_prior_infection}_REDACTED.html"),
      path = here::here("output", "report", "tables")
    )
  
}

################################################################################
# create a dataset for each analysis

# data_main for cancer vs noncancer
data_main <- data_tables %>%
  mutate(
    subgroup = if_else(cancer_subgroup == "noncancer", "noncancer", "cancer")
  ) 

# data_type for haem vs solid
data_type <- data_tables %>%
  mutate(
    subgroup = if_else(
      cancer_subgroup %in% c("cancer_haem", "cancer_solid"),
      str_remove(cancer_subgroup, "cancer_"),
      NA_character_
    )
  ) %>%
  filter(!is.na(subgroup)) 

# data_age for cancer vs non cancer in age groups
data_age <- data_tables %>%
  mutate(
    subgroup = str_c(if_else(cancer_subgroup == "noncancer", "noncancer", "cancer"), age_subgroup, sep = "; ")
  ) 

################################################################################
# apply the function to each of the analysis datasets
for (x in c(TRUE, FALSE)) {
  make_table1(data = data_main, include_prior_infection = x)
  make_table1(data = data_type, include_prior_infection = x)
  make_table1(data = data_age, include_prior_infection = x)
}
