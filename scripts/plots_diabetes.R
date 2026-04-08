library(tidyr)
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)
library(CI.RCT.Sim)

names <- c(
  "CS", # Core scenario
  "NC", # No correlation
  "RTE", # Reduced treatment effect
  "MR", # More rescue
  "PV", # Positivity violation
  "ROT", # Rescue on top
  "NM", # No missingness
  "HM", # Higher missingness
  "CS_null", # Core scenario - Null hypothesis
  "NC_null", # No correlation - Null hypothesis
  "RTE_null", # Reduced treatment effect - Null hypothesis
  "MR_null", # More rescue - Null hypothesis
  "PV_null", # Positivity violation - Null hypothesis
  "ROT_null", # Rescue on top - Null hypothesis
  "NM_null", # No missingness - Null hypothesis
  "HM_null" # Higher missingness - Null hypothesis
)


# Are there any scenarios that did not work?
# First check if it is a column name at all,
# second, filter out lists from the assumption data frame, since is.na doesn't work on lists
if ("FATAL_TERMINATION" %in% colnames(results)) {
  lists <- !sapply(results, is.list)
  results[, lists] |>
    filter(!is.na(FATAL_TERMINATION))
}

results_long <- results |>
  (\(dat){
    tmp <- attr(dat, "design_names")
    tmp$sim <- tmp$sim[!(tmp$sim == "FATAL_TERMINATION")]
    attr(dat, "design_names") <- tmp
    dat$FATAL_TERMINATION <- NULL
    dat
  })() |>
  mutate(
    scenario_nr = 1:n(),
    scenario_name = names
  ) |>
  results_pivot_longer()


# Plots for Estimators ----------------------------------------------------

results_long |>
  filter(str_detect(method, "\\.est")) |>
  ggplot(aes(x = scenario_nr, colour = method)) +
  aes(y = bias) +
  # geom_line() +
  geom_point() +
  geom_hline(yintercept = 0)

results_long |>
  filter(str_detect(method, "\\.est")) |>
  ggplot(aes(x = scenario_nr, colour = method)) +
  aes(y = coverage) +
  # geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.95) +
  scale_y_continuous(limits = c(0, 1))

results_long |>
  filter(str_detect(method, "\\.est")) |>
  ggplot(aes(x = scenario_nr, colour = method)) +
  aes(y = width) +
  # geom_line() +
  geom_point() +
  coord_cartesian(ylim = c(0, 4))


# Plots for Tests ---------------------------------------------------------

results_long |>
  filter(str_detect(method, "\\.test")) |>
  ggplot(aes(x = scenario_nr, colour = method)) +
  aes(y = rejection_0.05) +
  # geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.05) +
  scale_y_continuous(limits = c(0, 1))

results_long |>
  filter(str_detect(method, "\\.est")) |>
  ggplot(aes(x = scenario_nr, colour = method)) +
  aes(y = 1 - null_lower) +
  # geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.05) +
  scale_y_continuous(limits = c(0, 1))
