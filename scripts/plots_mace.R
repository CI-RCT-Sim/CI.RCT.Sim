library(tidyr)
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)
library(CI.RCT.Sim)

# Are there any scenarios that did not work?
results |>
  filter(!is.na(FATAL_TERMINATION))

results_long <- results |>
  (\(dat){
    tmp <- attr(dat, "design_names")
    tmp$sim <- tmp$sim[!(tmp$sim=="FATAL_TERMINATION")]
    attr(dat, "design_names") <- tmp
    dat$FATAL_TERMINATION <- NULL
    dat
  })() |>
  mutate(
    scenario_nr = 1:n()
  ) |>
  results_pivot_longer()


### List of abbreviations :
# BH : Baseline Hazard
# TE : Treatment Effect
# CE : Covariate Effect
# Disc : discontinuation
# Withd : Withdrawal
###
results_long <- results_long |> mutate(scenario_nr = case_when(
  scenario_nr==1 ~ "Baseline",
  scenario_nr==2 ~ "Higher MACE BH",
  scenario_nr==3 ~ "No TE High Sample Size",
  scenario_nr==4 ~ "Higher TE",
  scenario_nr==5 ~ "No TE buffer",
  scenario_nr==6 ~ "Lower TE buffer",
  scenario_nr==7 ~ "No CE MACE",
  scenario_nr==8 ~ "Confounding CE MACE",
  
  scenario_nr==9 ~ "Lower Disc BH",
  scenario_nr==10 ~ "No CE Disc",
  scenario_nr==11 ~ "Confounding CE Disc",
  
  scenario_nr==12 ~ "Higher Withd",
  scenario_nr==13 ~ "No CE Withd",
  scenario_nr==14 ~ "Confounding CE Withd",
  
  scenario_nr==15 ~ "Lower assumed window",
  scenario_nr==16 ~ "Higher Assumed window",
  
  scenario_nr==17 ~"No TE Low Sample Size"
))

# Plots for Estimators ----------------------------------------------------

results_long |>
  filter(str_detect(method, "\\.est")) |>
  ggplot(aes(x=scenario_nr, colour=method)) +
  aes(y=bias) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept=0)

results_long |>
  filter(str_detect(method, "\\.HR")) |>
  ggplot(aes(x=scenario_nr, colour=method)) +
  aes(y=bias) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept=0)

results_long |>
  filter(str_detect(method, "\\.HR")) |>
  ggplot(aes(x=scenario_nr, colour=method)) +
  aes(y=coverage) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept=0.95) +
  scale_y_continuous(limits=c(0,1))

results_long |>
  filter(str_detect(method, "\\.HR")) |>
  ggplot(aes(x=scenario_nr, colour=method)) +
  aes(y=width) +
  geom_line() +
  geom_point()


# Plots for Tests ---------------------------------------------------------

results_long |>
  filter(str_detect(method, "\\.test")) |>
  ggplot(aes(x=scenario_nr, colour=method)) +
  aes(y=rejection_0.05) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept=0.05) +
  scale_y_continuous(limits=c(0,1))

results_long |>
  filter(str_detect(method, "\\.HR")) |>
  ggplot(aes(x=scenario_nr, colour=method)) +
  aes(y=null_upper) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept=0.05) +
  scale_y_continuous(limits=c(0,1))

