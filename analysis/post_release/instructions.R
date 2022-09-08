# define release folder
release_folder <- here::here("release20220622")

### RESTART R SESSION IN BETWEEN RUNNING EACH SCRIPT

### Flow chart
# process data for flow chart
source(here::here("analysis", "post_release", "flow.R")) 

### Table 1
# process data for summary table in manuscript and supplementary material
source(here::here("analysis", "post_release", "table1_process"))
# render table for manuscript
relase_folder <- here::here("release20220622")
rmarkdown::render(
  here::here("analysis","post_release", "table1_process.Rmd"),
  knit_root_dir = release_folder,
    output_file = here::here(release_folder, "table1_process.docx"))

### Subsequent vaccination
# plot cumulative incidence of subsequent vaccination
source(here::here("analysis", "subsequent_vax", "plot_cumulative_incidence.R"))

### Metaregression
source(here::here("analysis", "post_release", "waning_metareg.R"))

### Results plots
# absolute risk
# source(here::here("analysis", "post_release", "plot_ar.R"))

# all hazard ratio plots
source(here::here("analysis", "post_release", "plot_cox_all.R"))
