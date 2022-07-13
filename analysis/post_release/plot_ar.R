library(tidyverse)
library(glue)
library(RColorBrewer)

release_folder <- here::here("release20220707")

ar_arm_weight <- readr::read_csv(file.path(release_folder, "ar_arm_weight.csv"))
ar_agesexonly <- readr::read_csv(file.path(release_folder, "ar_agesexonly.csv"))
ar_arm_group <- readr::read_csv(file.path(release_folder, "ar_arm_group.csv"))


outcomes_long <- c(
  "Positive SARS-CoV-2 test",
  "COVID-19 hospitalisation",
  "COVID-19 death"
  )

outcomes <- c("postest", "covidadmitted", "coviddeath")
names(outcomes_long) <- outcomes

cohorts <- c("cancer", "noncancer")
cohorts_long <- c("Cancer", "General population")

k_vals = 1:6#c(1,3,6)

data_clean <- ar_arm_group %>%
  filter(k %in% k_vals, outcome %in% outcomes) %>%
  select(-ar_crude) %>%
  mutate(across(ar_weighted, round, 5)) %>%
  mutate(across(outcome,
                factor,
                levels = outcomes,
                labels = str_wrap(outcomes_long[outcomes], 10))
  ) %>%
  mutate(across(cohort, 
                factor,
                levels = cohorts,
                labels = cohorts_long)
  ) %>%
  mutate(across(arm, 
                factor,
                levels = c("vax", "unvax"),
                labels = c("Vaccinated", "Unvaccinated"))
  ) %>%
  arrange(cohort, arm, outcome) 

palette_subgroups <- brewer.pal(n=4, name="Set2")[1:2]
names(palette_subgroups) <- cohorts_long

data_clean %>%
  ggplot(aes(x = k, y = ar_weighted, colour = cohort)) +
  geom_point(alpha = 0.8) +
  facet_grid(
    outcome ~ arm, 
    scales = "free_y",
    switch = "y"
    ) +
  labs(x = "Comparison Period", y = "Standardised absolute risk") +
  scale_x_continuous(
    breaks = 1:6
  ) +
  scale_colour_manual(
    values = palette_subgroups
  ) +
  theme_bw() +
  theme(
    
    panel.border = element_blank(),
    panel.spacing = unit(0.8, "lines"),
    # axis.line.y = element_line(colour = "black"),
    
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
    
    legend.position = "bottom"
  )


ggsave(
  filename = file.path(release_folder, glue("ar_plot.png")),
  width=16, height=14, units="cm"
  )

