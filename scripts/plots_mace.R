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
  scenario_nr==1 ~ "  Baseline",
  scenario_nr==2 ~ " MACE - Higher BH",
  scenario_nr==3 ~ "Null - High Sample Size",
  scenario_nr==4 ~ " TE - Higher",
  scenario_nr==5 ~ " TE - No TE buffer",
  scenario_nr==6 ~ " TE - Lower TE buffer",
  scenario_nr==7 ~ " MACE - No CE",
  scenario_nr==8 ~ " MACE - Confounding CE",
  
  scenario_nr==9 ~ " Disc - Lower BH",
  scenario_nr==10 ~ " Disc - No CE",
  scenario_nr==11 ~ " Disc - Confounding CE",
  
  scenario_nr==12 ~ " Withd - Higher BH",
  scenario_nr==13 ~ " Withd - No CE",
  scenario_nr==14 ~ " Withd - Confounding CE",
  
  scenario_nr==15 ~ " Assumed window - Lower",
  scenario_nr==16 ~ " Assumed window - Higher",
  scenario_nr==17 ~"Null - Confounding CE MACE",
  scenario_nr==18 ~"Null - Confounding CE Disc",
  scenario_nr==19 ~"Null - Confounding CE Withd",
  
  scenario_nr==20 ~"Null - Low Sample Size"
  
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

