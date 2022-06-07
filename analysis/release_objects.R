
# # # # # # # # # # # # # # # # # # # # #
# Purpose: To gather level 4 files ("moderately sensitive") place in a single directory for easy review and release
# # # # # # # # # # # # # # # # # # # # #

## Import libraries ----
library(tidyverse)
library(here)
library(glue)
library(survival)

release_objects <- "release_objects"
output_dir <- here("output", release_objects)
fs::dir_create(output_dir)
fs::dir_create(here(output_dir, "svp_plots"))

################################################################################
svp_plots <- list.files(path = here("output", "second_vax_period", "images"), 
                        pattern = "plot_by_region.+",
                        all.files = FALSE,
                        full.names = FALSE, recursive = FALSE,
                        ignore.case = FALSE, include.dirs = FALSE)

# SVP plots
for (i in svp_plots) {
  
  fs::file_copy(here("output", "second_vax_period", "images", i), 
                here("output", release_objects, "svp_plots", i), overwrite = TRUE)
  
}

# min and max follow-up dates for plots
fs::file_copy(here("output", "lib", "data_min_max_fu.csv"), 
              here("output", release_objects, "data_min_max_fu.csv"), overwrite = TRUE)

# survtable for subsequent vaccination
fs::file_copy(here("output", "subsequent_vax", "tables", "survtable_redacted.csv"), 
              here("output", release_objects, "survtable_redacted.csv"), overwrite = TRUE)

for (i in 1:4) {
  # table1
  fs::file_copy(here("output", "report", "tables", glue("table1_{i}_REDACTED.csv")), 
                here("output", release_objects, glue("table1_{i}_REDACTED.csv")), overwrite = TRUE)
  
}

# eligibility counts
fs::file_copy(here("output", "tables", "eligibility_count_all.csv"), 
              here("output", release_objects, "eligibility_count_all.csv"), overwrite = TRUE)

# event_counts
fs::file_copy(here("output", "models_cox", "data", "event_counts_all.csv"), 
              here("output", release_objects, "event_counts_all.csv"), overwrite = TRUE)

# estimates_all
fs::file_copy(here("output", "models_cox", "data", "estimates_all.csv"), 
              here("output", release_objects, "estimates_all.csv"), overwrite = TRUE)

################################################################################
## create text for output review issue ----
fs::dir_ls(here("output", "release_objects"), type="file", recurse=TRUE) %>%
  map_chr(~str_remove(., fixed(here()))) %>%
  map_chr(~paste0("- [ ] ", str_remove(.,fixed("/")))) %>%
  paste(collapse="\n") %>%
  writeLines(here("output", "files_for_release.txt"))
