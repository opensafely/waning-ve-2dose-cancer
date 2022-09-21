### RESTART R SESSION IN BETWEEN RUNNING EACH SCRIPT

### Flow chart
# process data for flow chart
source(here::here("analysis", "post_release", "flow.R")) 

### Subsequent vaccination
# plot cumulative incidence of subsequent vaccination
source(here::here("analysis", "subsequent_vax", "plot_cumulative_incidence.R"))

### Metaregression
source(here::here("analysis", "post_release", "waning_metareg.R"))

### Tables for manuscript / supplement

## Summarise characteristics
# process data for summary table in manuscript and supplementary material
source(here::here("analysis", "post_release", "table1_process"))
# render table for manuscript
relase_folder <- here::here("release20220622")
rmarkdown::render(
  here::here("analysis","post_release", "table1_process.Rmd"),
  knit_root_dir = release_folder,
  output_file = here::here(release_folder, "table1_process.docx"))

## tabulate hazard ratios

## tabulate ratios of hazard ratios 
source(here::here("analysis", "post_release", "rhr_table.R"))

### Figures for mauscript / supplement
# absolute risk
# source(here::here("analysis", "post_release", "plot_ar.R"))

# all hazard ratio plots
source(here::here("analysis", "post_release", "plot_cox_all.R"))
