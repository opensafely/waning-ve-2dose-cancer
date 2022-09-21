library(tidyverse)

release_folder <- here::here("release20220622")

estimates_all <- readr::read_csv(
  here::here(release_folder, "estimates_all.csv")
)

################################################################################
# read study parameters
study_parameters <- readr::read_rds(
  here::here("analysis", "lib", "study_parameters.rds"))

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))

# read outcomes
outcomes <- readr::read_rds(
  here::here("analysis", "lib", "outcomes.rds"))
outcomes_order <- c(2:5)

# scale for x-axis
K <- study_parameters$K
ends <- seq(2, (K+1)*4, 4)
starts <- ends + 1
weeks_since_2nd_vax <- str_c(starts[-(K+1)], ends[-1], sep = "-")

################################################################################

table_data <- estimates_all %>%
  filter(
    variable == "k",
    !reference_row,
    prior,
    outcome != "anytest",
    model == "max_adjusted"
  ) %>%
  mutate(across(c(estimate, starts_with("conf")), ~round(exp(.x), 2))) %>%
  transmute(
    subgroup = factor(subgroup, labels = subgroups),
    Brand = factor(
      comparison, 
      levels = c("BNT162b2", "ChAdOx1", "both"), 
      labels = c("BNT162b2", "ChAdOx1", "Combined")
      ),
    outcome,
    `Weeks since second dose` = factor(period, levels = 1:K, labels = weeks_since_2nd_vax), 
    value = glue::glue("{format(estimate, nsmall=2)} ({format(conf.low, nsmall=2)}, {format(conf.high, nsmall=2)})")
  ) %>%
  pivot_wider(
    names_from = outcome,
    values_from = value
  ) %>%
  mutate(across(any_of(unname(outcomes)),
                ~if_else(is.na(.x), "-", as.character(.x)))) %>%
  arrange(subgroup, Brand, `Weeks since second dose`) %>%
  select(subgroup, Brand, `Weeks since second dose`, outcomes[outcomes_order])
  
  
table_data_list <- table_data %>% 
  group_split(subgroup)

doc <- officer::read_docx() 

for (i in seq_along(table_data_list)) {
  
  ftab <- table_data_list[[i]] %>%
    select(-subgroup) %>%
    flextable::flextable() %>%
    flextable::merge_v(j=1, part = "body") %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::theme_booktabs() 
  
  doc <- doc %>%
    officer::body_add_caption(officer::block_caption(subgroups[i], style = "Normal")) %>%
    flextable::body_add_flextable(value = ftab, split = FALSE) %>%
    officer::body_add_break()
  
}

doc <- doc %>%
  print(target = here::here(release_folder, glue::glue("hr_tables.docx")))

