################################################################################
# This script:
# - combines estimates from all models into csv for release

################################################################################
library(tidyverse)
library(RColorBrewer)
library(lubridate)
library(glue)

## import command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  metareg=TRUE
} else{
  metareg=FALSE
}

################################################################################
release_folder <- here::here("release20220622")
if (!exists("release_folder")) {
  release_folder <- here::here("output", "release_objects")
  fs::dir_create(file.path(release_folder))
}
fs::dir_create(file.path(release_folder, "checking"))


################################################################################
# read study parameters
study_parameters <- readr::read_rds(
  here::here("analysis", "lib", "study_parameters.rds"))

# read outcomes
outcomes <- readr::read_rds(
  here::here("analysis", "lib", "outcomes.rds")
)
outcomes <- outcomes[outcomes!="covidemergency"]
outcomes_order <- c(3,4,2,5,1)
outcomes_long <- names(outcomes)
outcomes_long[outcomes=="covidadmitted"] <- "COVID-19 hospitalisation"
names(outcomes) <- outcomes_long
rm(outcomes_long)

# read subgroups
subgroups <- readr::read_rds(
  here::here("analysis", "lib", "subgroups.rds"))
subgroup_labels <- seq_along(subgroups)

# define comparisons
comparisons_old <- c("BNT162b2", "ChAdOx1", "both")
comparisons_new <- c("BNT162b2", "ChAdOx1", "Combined")

################################################################################
# read estimates data
estimates_all <- readr::read_csv(file.path(release_folder, "estimates_all.csv")) 

################################################################################
# gg plot pallete
gg_color_hue <- function(n, transparency = 1) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha = transparency)[1:n]
}

# function for adding transparency
addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

# function for estimated effectiveness axis
formatpercent100 <- function(x,accuracy){
  formatx <- scales::label_percent(accuracy)(x)
  
  if_else(
    formatx==scales::label_percent(accuracy)(1),
    paste0(">",scales::label_percent(1)((100-accuracy)/100)),
    formatx
  )
}

# scale for x-axis
K <- study_parameters$K
ends <- seq(2, (K+1)*4, 4)
starts <- ends + 1
weeks_since_2nd_vax <- str_c(starts[-(K+1)], ends[-1], sep = "-")

################################################################################
# page dimensions
page_height <- 26
page_width <- 17

# axis labels 
y_lab <- "Hazard Ratio (HR)"
y_lab_2 <- "Estimated vaccine effectiveness = 100*(1-HR)\n "
y_lab_adj <- "adjusted Hazard Ratio (aHR)"
y_lab_adj_2 <- "Estimated vaccine effectiveness = 100*(1-aHR)\n "
x_lab <- "Weeks since second dose"

# legend options
legend_width <- 15
subgroup_plot_labels <- subgroups

# derive data
plot_data <- estimates_all %>%
  filter(
    !reference_row,
    !is.na(estimate),
    variable %in% "k",
    # keep only outcomes of interest
    outcome %in% outcomes
  ) %>%
  mutate(k=as.integer(label)) %>%
  mutate(across(c(estimate, conf.low, conf.high), exp)) %>%
  mutate(across(subgroup,
                factor,
                levels = subgroup_labels,
                labels = subgroups
  )) 

# add metaregression data if available
if (metareg) {
  
  # read metareg data
  metareg_results_k <- readr::read_rds(
    file.path(release_folder, "metareg_results_k.rds")) %>%
    mutate(across(model,
                  ~case_when(
                    .x==1 ~ "unadjusted",
                    .x==2 ~ "part_adjusted",
                    .x==3 ~ "max_adjusted",
                    TRUE ~ NA_character_
                  ))) %>%
    select(model, subgroup, comparison, outcome, prior, k, starts_with("line")) 
  
  plot_data <- plot_data %>%
    left_join(
      metareg_results_k, 
      by = c("model", "subgroup", "comparison", "outcome", "prior", "k")
    ) %>%
    mutate(line_group = str_c(model, subgroup, comparison, outcome, prior, sep = "; ")) %>%
    # only plot line within range of estimates
    mutate(k_nonmiss = if_else(!is.na(estimate), k, NA_integer_)) %>%
    group_by(line_group) %>%
    mutate(keep = 0 < sum(!is.na(k_nonmiss))) %>%
    filter(keep) %>%
    mutate(
      min_k = min(k_nonmiss, na.rm = TRUE),
      max_k = max(k_nonmiss, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(across(line,
                  ~ if_else(min_k <= k & k <= max_k,
                            .x,
                            NA_real_)))
}

# add factor labels
plot_data <- plot_data %>%
  mutate(k_labelled = k) %>%
  mutate(across(k_labelled, 
                factor, 
                levels = 1:K,
                labels = weeks_since_2nd_vax)) %>%
  mutate(across(model,
                factor,
                levels = c("unadjusted", "part_adjusted", "max_adjusted"),
                labels = sapply(c("Stratfied Cox model, no further adjustment",
                                  "Stratfied Cox model, adjustment for demographic variables",
                                  "Stratfied Cox model, adjustment for demographic and clinical variables"),
                                str_wrap, width=100))) %>%
  mutate(outcome_unlabelled = outcome) %>%
  mutate(across(outcome,
                factor,
                levels = unname(outcomes[outcomes_order]),
                labels = str_wrap(names(outcomes[outcomes_order]), 10)
  )) %>%
  mutate(across(subgroup,
                factor,
                levels = subgroups,
                labels = str_wrap(subgroup_plot_labels, 100)
  )) %>%
  mutate(across(comparison,
                factor,
                levels = comparisons_old,
                labels = comparisons_new
  )) 

#################################################################################
# spacing of points on plot
position_dodge_val <- 0.6

# shape of points
shapes_subgroups <- c(21,22,23,24)
names(shapes_subgroups) <- subgroups[1:4]

linetypes_subgroups <- c("dotted", "longdash", "dashed", "dotdash")
names(linetypes_subgroups) <- subgroups[1:4]

# colour palettes
palette_subgroups <- brewer.pal(n=4, name="Set2")
palette_unadj <- addalpha(palette_subgroups, alpha=0.2)
palette_part <- addalpha(palette_subgroups, alpha=0.6)
palette_max <- palette_subgroups
names_models <- levels(plot_data$model)
names(palette_subgroups) <- subgroups[1:4]
names(palette_unadj) <- rep(names_models[1],4)
names(palette_part) <- rep(names_models[2],4)
names(palette_max) <- rep(names_models[3],4)

# breaks and lims for y-axes
primary_vax_y1 <- list(breaks = c(0.02, 0.05, 0.2, 0.5, 1, 2), 
                       limits = c(0.02, 2))
primary_vax_y2 <- list(breaks = c(0,0.5,0.8, 0.95, 0.98))
primary_brand_y1 <- list(breaks = c(0.2, 0.5, 1, 2, 5), 
                         limits = c(0.2, 5))
anytest_y1 <- list(breaks = c(0.5, 1, 2, 5), 
                   limits = c(0.5, 5))

# plot titles
plot_titles <- str_replace(subgroups, "_", ", ")
plot_titles <- str_replace(plot_titles, "noncancer", "General population cohort")  
plot_titles <- str_replace(plot_titles, "cancer", "Cancer cohort")
plot_titles <- str_replace(plot_titles, "cancer", "Cancer cohort")
plot_titles <- str_replace(plot_titles, "haem", "Haematological malignancy cohort")
plot_titles <- str_replace(plot_titles, "solid", "Solid organ malignancy")
plot_titles[str_detect(plot_titles, "\\d")] <- str_c(plot_titles[str_detect(plot_titles, "\\d")], " years")

################################################################################
# for each subgroup, compare all three models
# 1 plot per subgroup, facet_rows=outcome, facet_cols=brand, shades=model
plot_models <- function(group, include_prior_infection) {
  
  group_colour <- group
  if (group %in% 5:6) {
    group_colour <- 1
  } else if (group %in% 7:8) {
    group_colour <- 2
  }
  
  subtitle_str <- plot_titles[group]
  if (include_prior_infection) {
    subtitle_str <- str_c(subtitle_str, ", prior infection included")
  } else {
    subtitle_str <- str_c(subtitle_str, ", prior infection excluded")
  }
  
  palette_model <- c(
    palette_unadj[group_colour], 
    palette_part[group_colour], 
    palette_max[group_colour]
    )
  shape_subgroup <- shapes_subgroups[group_colour]
  linetype_subgroup <- linetypes_subgroups[group_colour]
  
  p <- plot_data %>%
    filter(
      subgroup == subgroups[group],
      prior == include_prior_infection,
      outcome_unlabelled != "anytest" 
    ) %>%
    droplevels() %>%
    ggplot(
      aes(
        x = k_labelled,
        colour = model,
        fill = model
      )
    ) +
    geom_hline(aes(yintercept=1), colour='grey') +
    geom_linerange(
      aes(ymin = conf.low, ymax = conf.high),
      position = position_dodge(width = position_dodge_val)
    ) +
    geom_point(
      aes(
        y = estimate,
        ),
      shape = shape_subgroup,
      position = position_dodge(width = position_dodge_val)
    ) +
    facet_grid(
      outcome ~ comparison, 
      switch = "y", 
      # scales = "free", 
      space = "free_x"
      ) +
    scale_y_log10(
      name = y_lab_adj,
      breaks = primary_vax_y1[["breaks"]],
      limits = primary_vax_y1[["limits"]],
      oob = scales::oob_keep,
      sec.axis = sec_axis(
        ~(1-.),
        name=y_lab_adj_2,
        breaks = primary_vax_y2[["breaks"]],
        labels = function(x){formatpercent100(x, 1)}
      )
    ) +
    labs(
      x = x_lab,
      subtitle = subtitle_str
    ) +
    scale_color_manual(values = palette_model, name = NULL) +
    scale_fill_manual(values = palette_model, name = NULL) +
    scale_shape_manual(values = shape_subgroup, name = NULL, guide = "none") +
    guides(colour = guide_legend(
      title = NULL,
      override.aes = list(shape = shape_subgroup, linetype = linetype_subgroup)
    )) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      axis.line.y = element_line(colour = "black"),
      
      axis.text = element_text(size=10),
      
      axis.title.x = element_text(size=10, margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(size=10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
      axis.text.x = element_text(size=8),
      axis.text.y = element_text(size=8),
      
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0),
      strip.text = element_text(size=8),
      
      panel.spacing = unit(0.8, "lines"),
      
      plot.title = element_text(hjust = 0),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.caption = element_text(hjust = 0, face= "italic"),
      
      legend.position = "bottom",
      legend.direction = "vertical",
      # big margins to cover up grid lines
      # legend.margin = margin(t = 30, r = 20, b = 30, l = 10),
      legend.key.width = unit(2, 'cm'),
      # legend.position = "bottom",
      legend.text = element_text(size=10)
    ) 
  
  ggsave(p,
         filename = file.path(release_folder, "checking", glue("hr_models_{group}_{include_prior_infection}.png")),
         width=page_height, height=page_width, units="cm")
    
}

for (x in 1:8) {
  for (y in c(TRUE, FALSE)) {
    try(plot_models(group=x, include_prior_infection = y))
  }
}

################################################################################
# for a given model, compare across subgroups, prior infection, age
plot_subgroups <- function(group, include_prior_infection, model) {
  

  # create subgroup_4 and age_group from subgroup
  subgroup_4_levels <- levels(plot_data$subgroup)[1:4]
  tmp_data <- plot_data %>%
    mutate(
      subgroup_4 = str_remove(as.character(subgroup), "_.+"),
      age_group = str_remove(as.character(subgroup), ".+_"),
        ) %>%
    mutate(across(subgroup_4, factor, levels = subgroup_4_levels)) %>%
    mutate(across(age_group, factor, levels = c("18-69", "70+")))
  
  # # filter subgroups and prior infection
  # filter_str <- sapply(
  #   seq_along(group),
  #   function(x)
  #     glue("(subgroup == subgroups[{group[x]}] & prior == include_prior_infection[{x}])")
  # )
  # filter_str <- str_c(filter_str, collapse = " | ")
  # filter_expr <- rlang::parse_expr(filter_str)
  # 
  
  # filter data according to arguments
  models <- levels(plot_data$model)
  m <- model
  tmp_data <- tmp_data %>%
    filter(
      subgroup %in% subgroups[group],
      prior %in% include_prior_infection,
      # !! filter_expr,
      model %in% models[m],
      !(outcome_unlabelled %in% c("anytest")) 
    ) %>%
    droplevels()
  
  subgroup_levels_current <- levels(tmp_data$subgroup_4)
  
  # indexes for colour palette
  group_index <- which(subgroup_4_levels %in% subgroup_levels_current)
  
  # indexes for fill
  age_group_levels <- levels(tmp_data$age_group)
  prior_levels <- unique(tmp_data$prior)
  # fill by prior infection unless more than one age group
  fill_by <- "prior"
  fill_var_levels <- as.character(sort(prior_levels))
  # white_level <- "FALSE" # when prior infection removed, white fill
  if (length(age_group_levels) > 1) {
    # if more than one age group, fill by age group
    fill_by <- "age_group"
    # white_level <- "18-69" # when age 18-69, white fill
    fill_var_levels <- age_group_levels
  } 
  
  # define colour palette
  palette_colour <- palette_subgroups[group_index]
  # define shapes and linetypes
  palette_shape <- shapes_subgroups[group_index]
  palette_linetype <- linetypes_subgroups[group_index]
  
  # define levels
  fill_var_levels <- unlist(lapply(
    seq_along(fill_var_levels), 
    function(x)
      str_c(subgroup_levels_current, fill_var_levels[x], sep = ", ")
  ))  
  fill_var_levels_clean <- str_replace(fill_var_levels, "TRUE", "prior infection included")
  fill_var_levels_clean <- str_replace(fill_var_levels_clean, "FALSE", "prior infection excluded")
  fill_var_levels_clean <- str_replace(fill_var_levels_clean, "18-69", "18-69 years")
  fill_var_levels_clean <- str_replace(fill_var_levels_clean, "70\\+", "70\\+ years")
  # derive variable  
  tmp_data <- tmp_data %>% 
    mutate(
      fill_var = factor(
        str_c(subgroup_4, !! sym(fill_by), sep = ", "),
        levels = fill_var_levels,
        labels = fill_var_levels_clean
      )
    )
  if (length(fill_var_levels) > 2) {
    # if using fill, add white to palette
    palette_fill <- c(rep("white",length(palette_colour)), palette_colour)
    names(palette_fill) <- fill_var_levels_clean
    legend_rows <- 2
  } else {
    # otherwise same as colour palette
    palette_fill <- palette_colour
    names(palette_fill) <- fill_var_levels_clean
    legend_rows <- 1
  }
  
  # create plot
  p <- tmp_data %>%
    ggplot(
      aes(
        x = k_labelled,
        colour = subgroup_4,
        fill = fill_var
      )
    ) +
    geom_hline(aes(yintercept=1), colour='grey') +
    geom_line(
      aes(y = line, 
          colour = subgroup_4, 
          linetype = subgroup_4,
          group = line_group
      )
    ) +
    geom_linerange(
      aes(ymin = conf.low, ymax = conf.high),
      position = position_dodge(width = position_dodge_val)
    ) +
    geom_point(
      aes(
        y = estimate,
        shape = subgroup_4
      ),
      position = position_dodge(width = position_dodge_val)
    ) +
    facet_grid(
      outcome ~ comparison, 
      switch = "y", 
      # scales = "free", 
      space = "free_x"
    ) +
    scale_y_log10(
      name = y_lab_adj,
      breaks = primary_vax_y1[["breaks"]],
      limits = primary_vax_y1[["limits"]],
      oob = scales::oob_keep,
      sec.axis = sec_axis(
        ~(1-.),
        name=y_lab_adj_2,
        breaks = primary_vax_y2[["breaks"]],
        labels = function(x){formatpercent100(x, 1)}
      )
    ) +
    labs(
      x = x_lab
    ) +
    scale_fill_manual(values = palette_fill, name = NULL) +
    scale_color_manual(values = palette_colour, name = NULL, guide="none") +
    scale_shape_manual(values = palette_shape, name = NULL, guide = "none") +
    scale_linetype_manual(values = palette_linetype, guide = "none") +
    guides(
      fill = guide_legend(
        title = NULL,
        nrow = legend_rows,
        override.aes = list(
          colour = rep(palette_colour,times=legend_rows),
          shape = rep(palette_shape,times=legend_rows),
          linetype = rep(palette_linetype,time=legend_rows),
          fill = palette_fill
        )
      )
    ) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      axis.line.y = element_line(colour = "black"),
      
      axis.text = element_text(size=10),
      
      axis.title.x = element_text(size=10, margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(size=10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
      axis.text.x = element_text(size=8),
      axis.text.y = element_text(size=8),
      
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0),
      strip.text = element_text(size=8),
      
      panel.spacing = unit(0.8, "lines"),
      
      plot.title = element_text(hjust = 0),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.caption = element_text(hjust = 0, face= "italic"),
      
      legend.position = "bottom",
      # legend.direction = "vertical",
      # big margins to cover up grid lines
      # legend.margin = margin(t = 30, r = 20, b = 30, l = 10),
      legend.key.width = unit(2, 'cm'),
      # legend.position = "bottom",
      legend.text = element_text(size=10)
    ) 
  
  filename_group <- str_c(group, collapse = "")
  filename_prior <- include_prior_infection
  if (length(include_prior_infection)>1) filename_prior <- "prior"
  filename <- glue("hr_subgroups_{filename_group}_{filename_prior}_{model}.png")
  
  ggsave(p,
         filename = file.path(release_folder, "checking", filename),
         width=page_height, height=page_width, units="cm")
  
}

# subgroups 1:2, no fill
## cancer vs noncancer, main comparison (Lee)
try(plot_subgroups(group=1:2, include_prior_infection=TRUE, model=1))
## haem vs solid, main comparison (Lee)
try(plot_subgroups(group=3:4, include_prior_infection=TRUE, model=2))

# subgroups 1:2, fill by prior
try(plot_subgroups(group=1:2, include_prior_infection=c(FALSE,TRUE), model=2))

# subgroups 1:2, fill by age i.e. subgroups 5:8)
try(plot_subgroups(group=5:8, include_prior_infection=TRUE, model=3))


##############################################################################
# 3-6 month estimates from Lee study
lee_data <- tribble(
  ~subgroup_4, ~comparison, ~prior, ~outcome_unlabelled, ~model, ~lee_period, ~estimate, ~conf.low, ~conf.high, 
  "noncancer", "Combined", TRUE, "postest", 1, "overall", 61.4, 61.4, 61.5,
  "noncancer", "Combined", TRUE, "postest", 1, "3-6 months", 69.8, 69.8, 69.9,
  "noncancer", "Combined", TRUE, "postest", 1, "0-8 weeks", 80.7, NA_real_, NA_real_,
  "noncancer", "Combined", TRUE, "postest", 1, "8-16 weeks", 64.9, NA_real_, NA_real_,
  "noncancer", "Combined", TRUE, "postest", 1, "16-24 weeks", 61.2, NA_real_, NA_real_,
  "noncancer", "Combined", TRUE, "postest", 1, "24-32 weeks", 69.6, NA_real_, NA_real_,
  "cancer", "Combined", TRUE, "postest", 1, "overall", 65.5, 65.1, 65.9,
  "cancer", "Combined", TRUE, "postest", 1, "3-6 months", 47.0, 46.3, 47.6,
  "cancer", "Combined", TRUE, "postest", 1, "0-8 weeks", 90.5, NA_real_, NA_real_,
  "cancer", "Combined", TRUE, "postest", 1, "8-16 weeks", 66.8, NA_real_, NA_real_,
  "cancer", "Combined", TRUE, "postest", 1, "16-24 weeks", 42.5, NA_real_, NA_real_,
  "cancer", "Combined", TRUE, "postest", 1, "24-32 weeks", 35.4, NA_real_, NA_real_,
  "cancer", "Combined", TRUE, "covidadmitted", 1, "overall", 84.5, 83.6, 85.4,
  "cancer", "Combined", TRUE, "covidadmitted", 1, "3-6 months", 74.6, 72.8, 76.3,
  "cancer", "Combined", TRUE, "coviddeath", 1, "overall", 93.5, 93.0, 94.0,
  "cancer", "Combined", TRUE, "coviddeath", 1, "3-6 months", 90.3, 89.3, 91.2
)

summary_plot <- local({
  
  # arguments
  group <- 1:2
  include_prior_infection <- TRUE
  model <- 1
  
  # create subgroup_4 and age_group from subgroup
  subgroup_4_levels <- levels(plot_data$subgroup)[1:4]
  tmp_data <- plot_data %>%
    mutate(
      subgroup_4 = str_remove(as.character(subgroup), "_.+"),
      age_group = str_remove(as.character(subgroup), ".+_"),
    ) %>%
    mutate(across(subgroup_4, factor, levels = subgroup_4_levels)) %>%
    mutate(across(age_group, factor, levels = c("18-69", "70+")))
  
  # filter data according to arguments
  models <- levels(plot_data$model)
  m <- model
  tmp_data <- tmp_data %>%
    filter(
      comparison == "Combined",
      subgroup %in% subgroups[group],
      prior %in% include_prior_infection,
      # !! filter_expr,
      model %in% models[m],
      !(outcome_unlabelled %in% c("anytest", "noncoviddeath")) 
    ) %>%
    select(-comparison) %>%
    bind_rows(
      lee_data %>%
        filter(
          lee_period %in% c("3-6 months", "overall")
        ) %>%
        mutate(
          xmin_val = if_else(
            lee_period == "3-6 months", 3.5, 1
          ),
          xmax_val = 6,
          ) %>%
        select(-model, -comparison) %>%
        mutate(across(subgroup_4,
                      factor,
                      levels = str_wrap(subgroup_plot_labels, 100)
        )) %>%
        mutate(outcome = outcome_unlabelled) %>%
        mutate(across(outcome,
                      factor,
                      levels = unname(outcomes[outcomes_order]),
                      labels = str_wrap(names(outcomes[outcomes_order]), 10)
        )) %>%
        # uncount(weights=3) %>%
        select(-lee_period) %>%
        mutate(k_labelled = 6) %>%
        mutate(across(k_labelled,
                      factor,
                      levels = 1:K,
                      labels = weeks_since_2nd_vax)) %>%
        mutate(across(c(estimate, conf.low, conf.high), ~ 1 - (.x/100))) %>%
        rename_with(.fn = ~ str_c("lee_", .x), .cols = c(estimate, conf.low, conf.high)) #%>%
        # mutate(line_group = str_c(subgroup, outcome_unlabelled, sep = "; "))
    ) %>%
    droplevels() 
  
  subgroup_levels_current <- levels(tmp_data$subgroup_4)
  
  # indexes for colour palette
  group_index <- which(subgroup_4_levels %in% subgroup_levels_current)
  
  # indexes for fill
  age_group_levels <- levels(tmp_data$age_group)
  prior_levels <- unique(tmp_data$prior)
  # fill by prior infection unless more than one age group
  fill_by <- "prior"
  fill_var_levels <- as.character(sort(prior_levels))
  # white_level <- "FALSE" # when prior infection removed, white fill
  if (length(age_group_levels) > 1) {
    # if more than one age group, fill by age group
    fill_by <- "age_group"
    # white_level <- "18-69" # when age 18-69, white fill
    fill_var_levels <- age_group_levels
  } 
  
  # define colour palette
  palette_colour <- palette_subgroups[group_index]
  # define shapes and linetypes
  palette_shape <- shapes_subgroups[group_index]
  palette_linetype <- linetypes_subgroups[group_index]
  
  # define levels
  fill_var_levels <- unlist(lapply(
    seq_along(fill_var_levels), 
    function(x)
      str_c(subgroup_levels_current, fill_var_levels[x], sep = ", ")
  ))  
  fill_var_levels_clean <- str_replace(fill_var_levels, "TRUE", "prior infection included")
  fill_var_levels_clean <- str_replace(fill_var_levels_clean, "FALSE", "prior infection excluded")
  fill_var_levels_clean <- str_replace(fill_var_levels_clean, "18-69", "18-69 years")
  fill_var_levels_clean <- str_replace(fill_var_levels_clean, "70\\+", "70\\+ years")
  # derive variable  
  tmp_data <- tmp_data %>% 
    mutate(
      fill_var = factor(
        str_c(subgroup_4, !! sym(fill_by), sep = ", "),
        levels = fill_var_levels,
        labels = fill_var_levels_clean
      )
    )
  if (length(fill_var_levels) > 2) {
    # if using fill, add white to palette
    palette_fill <- c(rep("white",length(palette_colour)), palette_colour)
    names(palette_fill) <- fill_var_levels_clean
    legend_rows <- 2
  } else {
    # otherwise same as colour palette
    palette_fill <- palette_colour
    names(palette_fill) <- fill_var_levels_clean
    legend_rows <- 1
  }
  
  # create plot
  p <- tmp_data %>%
    ggplot(
      aes(
        x = period,
        colour = subgroup_4,
        fill = fill_var
      )
    ) +
    geom_rect(
      # xmin = 3.7, xmax = 7,
      aes(
        xmin = xmin_val, xmax = xmax_val,
        ymax = lee_conf.low, ymin = lee_conf.high,
        fill = fill_var
      ),
      alpha = 0.2,
      # colour = NA,
      linetype = "solid"
    ) +
    geom_hline(aes(yintercept=1), colour='grey') +
    geom_line(
      aes(y = line, 
          colour = subgroup_4, 
          linetype = subgroup_4,
          group = line_group
      )
    ) +
    geom_linerange(
      aes(ymin = conf.low, ymax = conf.high),
      position = position_dodge(width = position_dodge_val)
    ) +
    geom_point(
      aes(
        y = estimate,
        shape = subgroup_4
      ),
      position = position_dodge(width = position_dodge_val)
    ) +
    facet_grid(
      outcome ~ ., 
      switch = "y", 
      # scales = "free", 
      space = "free_x"
    ) +
    scale_x_continuous(
      breaks = seq_along(weeks_since_2nd_vax),
      labels = weeks_since_2nd_vax,
      limits = c(1, length(weeks_since_2nd_vax))
    ) +
    scale_y_log10(
      name = y_lab_adj,
      breaks = primary_vax_y1[["breaks"]],
      limits = primary_vax_y1[["limits"]],
      oob = scales::oob_keep,
      sec.axis = sec_axis(
        ~(1-.),
        name=y_lab_adj_2,
        breaks = primary_vax_y2[["breaks"]],
        labels = function(x){formatpercent100(x, 1)}
      )
    ) +
    labs(
      x = x_lab
    ) +
    scale_fill_manual(values = palette_fill, name = NULL) +
    scale_color_manual(values = palette_colour, name = NULL, guide="none") +
    scale_shape_manual(values = palette_shape, name = NULL, guide = "none") +
    scale_linetype_manual(values = palette_linetype, guide = "none") +
    guides(
      fill = guide_legend(
        title = NULL,
        nrow = legend_rows,
        override.aes = list(
          colour = rep(palette_colour,times=legend_rows),
          shape = rep(palette_shape,times=legend_rows),
          linetype = rep(palette_linetype,time=legend_rows),
          fill = palette_fill
        )
      )
    ) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      axis.line.y = element_line(colour = "black"),
      
      axis.text = element_text(size=10),
      
      axis.title.x = element_text(size=10, margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(size=10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
      axis.text.x = element_text(size=8),
      axis.text.y = element_text(size=8),
      
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0),
      strip.text = element_text(size=8),
      
      panel.spacing = unit(0.8, "lines"),
      
      plot.title = element_text(hjust = 0),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.caption = element_text(hjust = 0, face= "italic"),
      
      legend.position = "bottom",
      # legend.direction = "vertical",
      # big margins to cover up grid lines
      # legend.margin = margin(t = 30, r = 20, b = 30, l = 10),
      legend.key.width = unit(2, 'cm'),
      # legend.position = "bottom",
      legend.text = element_text(size=10)
    ) 
  
  ggsave(p,
         filename = file.path(release_folder, "checking", "hr_subgroups_lee.png"),
         width=page_width, height=16, units="cm")
  
(boxes)})
