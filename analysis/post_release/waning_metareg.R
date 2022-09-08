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
data_metareg_0 <- estimates_all %>%
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

# expand to all combinations
data_metareg <- data_metareg_0 #%>%
  # expand(subgroup, comparison, outcome, model, prior, k) %>%
  # left_join(data_metareg_0) #%>%
  # # following lines because in a few instances seloghr=0, and this is breaking the metareg code 
  # mutate(seloghr = (conf.high - conf.low)/(2*1.96)) %>%
  # mutate(across(c(estimate, conf.low, conf.high),
  #               ~if_else(near(seloghr, 0), NA_real_, .x))) %>%
  # select(-seloghr)

################################################################################

data_metareg_split <- data_metareg %>%
  # filter(!is.na(estimate)) %>%
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
          # I2 describes the percentage of the variability in effect estimates that is due to heterogeneity rather than sampling error (chance).
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


# plot_data <- rma_res %>%
#   left_join(data_metareg %>%
#               filter(sex == "Both", ageband == "all") %>%
#               select(subgroup, comparison, outcome, k, estimate, conf.low, conf.high) %>%
#               mutate(across(c(estimate, conf.low, conf.high), exp)) %>%
#               mutate(across(k, ~.x+1)), 
#             by = c("subgroup", "comparison", "outcome")) %>%
#   mutate(
#     line = loghr1 + (k-1)*logrhr#, # use k-1 because intercept at k=-1
#     # CI for slope not valid as don't know covariance between intercept and slope
#     # line_lower = loghr1_lower + (k-1)*logrhr_lower,
#     # line_higher = loghr1_higher + (k-1)*logrhr_higher
#   ) %>%
#   mutate(across(starts_with("line"), exp)) %>%
#   mutate(across(outcome, factor, levels = c("covidadmitted", "coviddeath", "postest", "noncoviddeath", "anytest"))) %>%
#   mutate(across(subgroup, factor, subgroups)) %>%
#   rename_with(~str_c("bnt_", .x), ends_with("2")) %>%
#   mutate(
#     bnt_logrhr = logrhr, az_logrhr = bnt_logrhr,
#     az_tau2 = bnt_tau2, az_I2 = bnt_I2, az_H2 = bnt_H2, az_R2 = bnt_R2
#   ) %>%
#   mutate(across(starts_with("az_"), ~if_else(k==1 & comparison == "ChAdOx1", signif(.x,3), NA_real_))) %>%
#   mutate(across(starts_with("bnt_"), ~if_else(k==2 & comparison == "BNT162b2", signif(.x,3), NA_real_))) 
# 
# 
# 
# plot_data %>%
#   mutate(text_x = 4, text_y = 0.01) %>%
#   filter(comparison != "both") %>%
#   ggplot(aes(
#     x = k, 
#     colour = comparison, 
#     shape = comparison,
#     fill = comparison
#   )) +
#   geom_hline(aes(yintercept=1), colour='grey') +
#   geom_line(
#     aes(y = line, 
#         colour = comparison, 
#         linetype = comparison,
#         group = comparison), 
#     alpha = 0.6
#   ) +
#   geom_linerange(
#     aes(ymin = conf.low, ymax = conf.high)#,
#     # position = position_dodge(width = position_dodge_val)
#   ) +
#   geom_point(
#     aes(y = estimate)#,
#     # position = position_dodge(width = position_dodge_val)
#   ) +
#   geom_text(aes(x = text_x, y = text_y + 0.05, label = str_c("logrhr=", bnt_logrhr, "; tau2=", bnt_tau2, "; I2=", bnt_I2, "%; R2=", bnt_R2, "%"))) +
#   geom_text(aes(x = text_x, y = text_y, label = str_c("logrhr=", az_logrhr, "; tau2=", az_tau2, "; I2=", az_I2, "%; R2=", az_R2, "%"))) +
#   facet_grid(outcome ~ subgroup, switch = "y", space = "free_x") +
#   scale_y_log10(
#     # name = y_lab_adj,
#     # breaks = primary_vax_y1[["breaks"]],
#     # limits = primary_vax_y1[["limits"]],
#     oob = scales::oob_keep#,
#     # sec.axis = sec_axis(
#     #   ~(1-.),
#     #   # name=y_lab_adj_2,
#     #   # breaks = primary_vax_y2[["breaks"]],
#     #   # labels = function(x){formatpercent100(x, 1)}
#     # )
#   ) +
#   # labs(
#   #   x = x_lab
#   # ) +
#   # scale_fill_discrete(guide = "none") +
#   # scale_shape_manual(values = comparison_shapes[1:2], name = NULL) +
#   # scale_color_manual(values = palette_adj[1:2], name = NULL) +
#   # scale_linetype_manual(values = comparison_linetypes[1:2], name = NULL) +
#   # guides(shape = guide_legend(
#   #   title = NULL, 
#   #   override.aes = list(colour = palette_adj[1:2], fill = comparison_shapes[1:2])
#   # )) +
#   theme_bw()
# 
# ######
# 
# 
# rma_res %>%
#   ggplot(aes(x = comparison, y = tau2)) +
#   geom_point() +
#   facet_grid(subgroup ~ outcome) 
#   
# 
# rma_res %>%
#   select(-contains(c("log", "tau", "H2"))) %>%
#   pivot_longer(cols = c(I2, R2)) %>%
#   ggplot(aes(x = comparison, y = value, colour = name)) +
#   geom_point() +
#   facet_grid(subgroup ~ outcome) +
#   theme_bw()
# 
# 
# 
# # data_select <- data_metareg %>%
# #   filter(
# #     subgroup == "65+ years",
# #     comparison == "BNT162b2",
# #     sex == "Both",
# #     ageband == "all",
# #     outcome == "covidadmitted",
# #     model == "adjusted"
# #   )
# # 
# # # rma runs a random-effects meta-analysis, which extended to
# # # mixed-effects meta-regression models when moderators (mods arg) are added.
# # # options specified to match stata defaults
# # metafor_res <- metafor::rma(
# #   yi = estimate,
# #   sei = seloghr,
# #   mods = ~k,
# #   method = "REML", # restricted maximum likelihood estimator (Viechtbauer, 2005; Raudenbush, 2009)
# #   test = "knha", # Knapp-Hartung method
# #   data = data_select
# # )
# # 
# # meta::metareg
# # 
# # bind_rows(
# #   stata_select %>%
# #     mutate(package = "stata") %>%
# #     select(package, logrhr, selogrhr, loghr1, seloghr1),
# #   tibble(
# #     package = "metafor",
# #     logrhr = metafor_res$beta[2,1],
# #     selogrhr = metafor_res$se[2],
# #     loghr1 = metafor_res$beta[1,1],
# #     seloghr1 = metafor_res$se[1],
# #   )
# # )
# 
# 
# 
