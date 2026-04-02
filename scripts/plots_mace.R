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
