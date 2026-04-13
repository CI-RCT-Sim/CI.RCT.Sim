library(tidyr)
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)
library(CI.RCT.Sim)

# # Are there any scenarios that did not work?
# results |>
#   filter(!is.na(FATAL_TERMINATION))

results_long <- results |>
  (\(dat){
    tmp <- attr(dat, "design_names")
    tmp$sim <- tmp$sim[!(tmp$sim=="FATAL_TERMINATION")]
    attr(dat, "design_names") <- tmp
    dat$FATAL_TERMINATION <- NULL
    dat
  })() |>
  arrange(desc(beta_A2)) |>
  mutate(
    scenario_nr = factor(1:n())
  ) |>
  results_pivot_longer()


# Plots for Estimators ----------------------------------------------------

update_theme(
  legend.position="bottom"
)

results_long |>
  filter(str_detect(method, "\\.est")) |>
  ggplot(aes(x=scenario_nr, colour=method, group=interaction(method, beta_A2))) +
  aes(y=bias) +
  geom_line() +
  geom_point(size=2.5) +
  geom_hline(yintercept=0) +
  labs(
    title = "Bias for Estimators in the Vaccine Scenario",
    y = "Bias [VE]",
    x = "Scenario"
  )

results_long |>
  filter(str_detect(method, "\\.est")) |>
  ggplot(aes(x=scenario_nr, colour=method, group=interaction(method, beta_A2))) +
  aes(y=coverage) +
  geom_line() +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.95) +
  scale_y_continuous(limits=c(0,1)) +
  labs(
    title = "CI Coverage for Estimators in the Vaccine Scenario",
    y = "Coverage",
    x = "Scenario"
  )

# results_long |>
#   filter(str_detect(method, "\\.est")) |>
#   ggplot(aes(x=scenario_nr, colour=method, group=interaction(method, beta_A2))) +
#   aes(y=width) +
#   geom_line() +
#   geom_point() +
#   coord_cartesian(ylim=c(0, 4))


# Plots for Tests ---------------------------------------------------------

results_long |>
  mutate(Hnull = if_else(VE <= 0.3, "Null", "Alternative")) |>
  filter(str_detect(method, "\\.test")) |>
  ggplot(aes(x=scenario_nr, colour=method, group=interaction(method, beta_A2), shape=Hnull)) +
  aes(y=rejection_0.05) +
  geom_line() +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.05) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, by=0.1), minor_breaks = seq(0, 1, by=0.05)) +
  labs(
    title = "Rejection Rate for the Tests in the Vaccine Scenario",
    y = "Rejection [alpha=0.05]",
    x = "Scenario",

  )

results_long |>
  mutate(Hnull = if_else(VE <= 0.3, "Null", "Alternative")) |>
  filter(str_detect(method, "\\.est")) |>
  ggplot(aes(x=scenario_nr, colour=method, group=interaction(method, beta_A2), shape=Hnull)) +
  aes(y=1-null_lower) +
  geom_line() +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.025) +
  scale_y_continuous(limits=c(0,1)) +
  labs(
    title = "Rejection Rate for the CI Based Tests in the Vaccine Scenario",
    y = "Rejection [alpha=0.025]",
    x = "Scenario",
  )
