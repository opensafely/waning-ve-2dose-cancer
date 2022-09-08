library(tidyverse)
library(glue)
library(flextable)
library(officer)
library(magrittr)
library(kableExtra)

################################################################################
release_folder <- here::here("release20220622")

################################################################################
# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))

subgroups_long <- subgroups
subgroups_long <- str_replace(subgroups_long, "noncancer", "General population")
subgroups_long <- str_replace(subgroups_long, "cancer", "Cancer cohort")
subgroups_long <- str_replace(subgroups_long, "haem", "Haematological malignancy")
subgroups_long <- str_replace(subgroups_long, "solid", "Solid organ malignancy")
subgroups_long <- str_replace(subgroups_long, "_", " ")
subgroups_long <- str_remove(subgroups_long, "\\d+.+")

################################################################################
all_table1_files <- list.files(path = release_folder, 
                        pattern = "table1_.+.csv",
                        all.files = FALSE,
                        full.names = FALSE, recursive = FALSE,
                        ignore.case = FALSE, include.dirs = FALSE)

# only those with prior infection included
all_table1_files_TRUE <- all_table1_files[str_detect(all_table1_files, "TRUE")]

# prepare the table 1 data
table_out0_list <- lapply(
  all_table1_files_TRUE,
  function(x) 
    readr::read_csv(here::here(release_folder, x), show_col_types = FALSE) %>%
    mutate(across(starts_with(c("BNT162b2", "ChAdOx1", "Unvaccinated")),
                  ~ if_else(
                    is.na(.x),
                    "-",
                    .x))) %>%
    mutate(across(Variable, ~if_else(Characteristic == "N", "N", .x))) %>%
    mutate(across(Characteristic, ~if_else(.x == "N", "", .x))) 
) 

names(table_out0_list) <- str_remove_all(all_table1_files_TRUE, "table1_|_TRUE_REDACTED.csv")

table_out0_list$`18-69` <- table_out0_list$age %>%
  select(Variable, Characteristic, ends_with(c("5"))) %>%
  full_join(
    table_out0_list$age %>%
      select(Variable, Characteristic, ends_with(c("7")))
  )

table_out0_list$`70+` <- table_out0_list$age %>%
  select(Variable, Characteristic, ends_with(c("6"))) %>%
  full_join(
    table_out0_list$age %>%
      select(Variable, Characteristic, ends_with(c("8")))
  )

table_out0_list <- table_out0_list[!(names(table_out0_list)=="age")]

################################################################################

create_table <- function(table_out0) {
  
  title <- deparse(substitute(table_out0))
  title <- str_remove(title, "table_out0_list\\$")
  title <- str_remove_all(title, "`")
  
  ################################################
  # prepare the column names
  # column names for table
  table_names <- names(table_out0)
  # get subgroup labels
  subgroup_labels <- str_extract(str_extract(table_names, "_\\d"), "\\d")
  subgroup_labels <- unique(subgroup_labels[!is.na(subgroup_labels)])
  # remove subgroup label
  col_names <- str_remove(table_names, "_\\d")
  names(col_names) <- table_names
  # number of columns for each subgroup
  cols_subtype <- sapply(
    subgroup_labels, 
    function(x) sum(str_detect(table_names, glue("_{x}")))
  )
  names(cols_subtype) <- subgroups_long[as.integer(subgroup_labels)]
  # # reorder
  # cols_subtype <- cols_subtype[which(names(cols_subtype) %in% as.integer(subgroup_labels))]
  
  ################################################
  # create and save the version for the manuscript
  
  # clean characteristic column
  table_out1 <- table_out0 %>%
    mutate(across(Characteristic, ~str_remove(.x, "\\s\\w+\\sdeprived"))) %>%
    mutate(across(Characteristic, 
                  ~if_else(
                    Variable == "BMI",
                    str_remove_all(str_extract(.x, "\\(\\d.+\\)"), "\\(|\\)"),
                    .x
                  ))) %>%
    mutate(across(Characteristic, 
                  ~if_else(
                    Variable == "BMI" & is.na(.x),
                    "<30",
                    .x
                  ))) 
  
  # order characteristics
  table_out1_list <- table_out1 %>%
    group_split(Variable) %>%
    as.list()
  
  names_table_out1_list <- lapply(
    table_out1_list, 
    function(x) unique(x$Variable)
  )
  
  names(table_out1_list) <- names_table_out1_list
  
  table_out1_list$`JCVI group` <- table_out1_list$`JCVI group` %>%
    filter(!(Characteristic %in% c("99", "01"))) %>%
    arrange(Characteristic) 
  
  table_out1_list$BMI <- table_out1_list$BMI %>%
    arrange(Characteristic) 
  
  table_out1_ordered <- bind_rows(table_out1_list)
  
  variable_order_manuscript <- c(
    "N" = "N", 
    "Age" = "Age", 
    "Sex" = "Sex", 
    "JCVI group" = "JCVI group",
    "Region" = "Region",
    "IMD" = "IMD (1 is most deprived)", 
    "Ethnicity" = "Ethnicity",
    "BMI" = "BMI", 
    "Learning disability" = "Learning disability",
    "Serious mental illness" = "Serious mental illness",
    "Morbidity count" = "Morbidity count",
    "Number of SARS-CoV-2 tests between 2020-05-18 and min_elig_date" = "Number of SARS-CoV-2 tests", 
    "Prior SARS-CoV-2 infection" = "Prior SARS-CoV-2 infection",
    "Flu vaccine in previous 5 years" = "Flu vaccine",
    "Pregnancy" = "Pregnancy")
  
  col_names[1] <- "Characteristic"
  
  table1_data_manuscript <- tibble(
    variable_long = names(variable_order_manuscript),
    variable_short = unname(variable_order_manuscript)
  ) %>%
    left_join(table_out1_ordered, by = c("variable_long" = "Variable")) %>%
    select(-variable_long) %>%
    rename(Variable = variable_short)
  
  # create table1_manuscript.docx
  page_width_docx <- 26 #cm
  cell_padding <- 0 # this is a guess, not sure what default is
  col1_with <- 1.75
  col2_width <- 1.5
  coli_width <- (page_width_docx - cell_padding*ncol(table1_data_manuscript) - col1_with - col2_width)/(ncol(table1_data_manuscript) - 2)
  
  flextable1 <- table1_data_manuscript %>%
    flextable() %>%
    set_header_labels(
      values = as.list(col_names)
    ) %>%
    merge_v(j=1, part = "body") %>%
    merge_at(i=1,j=1:2, part="header") %>%
    add_header_row(
      values = c("", names(cols_subtype)),
      colwidths = c(2, unname(cols_subtype))
    ) %>%
    width(j=1, width = col1_with, unit = "cm") %>%
    width(j=2, width = col1_with, unit = "cm") %>%
    width(j=3:ncol(table1_data_manuscript), width = coli_width, unit = "cm") %>%
    fontsize(size = 8, part = "all") %>%
    theme_booktabs() %>%
    padding(
      padding.top = 1,
      padding.bottom = 1,
      part = "all"
    )
  
  return(flextable1)
  
  # doc <- read_docx() %>%
  #   body_add_flextable(value = flextable1, split = FALSE) %>%
  #   body_end_section_landscape() %>% # a landscape section is ending here
  #   print(target = here::here(release_folder, glue("table1_{title}.docx")))
  
}

create_table(table_out0_list$main)
create_table(table_out0_list$`18-69`)
create_table(table_out0_list$`70+`)
create_table(table_out0_list$type)
