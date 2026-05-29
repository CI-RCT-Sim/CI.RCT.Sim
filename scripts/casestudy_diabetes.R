source("scripts/plots_common.R")

cond <- diabetes_scenario(print = FALSE)[1,] |>
  diabetes_scenario_set_truevalues()
set.seed(201)
dat <- generate_diabetes(cond)

## Convert in to long format
dat_long <-
  dat |>
  pivot_longer(
    cols = matches("^(y|R)\\d+$"),
    names_to = c(".value", "visit"),
    names_pattern = "([yR])(\\d+)"
  ) |>
  mutate(
    visit = as.integer(visit),
    rescue = R,
    missing = is.na(y)
  )

## Descriptive tables with rescue medication use and missingness

desc_tbl <-
  dat_long |>
  group_by("Treatment" = trt,
           "Visit" = visit) |>
  summarise(
    n = n(),
    "Mean HbA1c" = mean(y, na.rm = TRUE),
    # sd_hba1c = sd(y, na.rm = TRUE),
    # rescue_n = sum(rescue == 1, na.rm = TRUE),
    "Rescue proportion" = 100 * mean(rescue == 1, na.rm = TRUE),
    # missing_n = sum(missing),
    "Missingness proportion" = 100 * mean(missing),
    .groups = "drop"
  )

desc_tbl

knitr::kable(
  desc_tbl,
  format = "latex",
  caption = "Descriptive table with visit specific rescue medication use and missingness",
  label = "diabetes_casestudy_descriptive"
) |>
  cat(file = "tables/diabetes_casestudy_descriptive.tex")


## Plots with trajectories of HbA1c values over time

trajectories <- ggplot(dat_long, aes(x = visit, y = y, group = id, colour = factor(trt))) +
  geom_line(alpha = 0.1) +
  # geom_point(alpha = 0.25) +
  stat_summary(aes(group = trt), fun = mean, geom = "line", linewidth = 1.2) +
  labs(x = "Visit", y = "HbA1c", colour = "Treatment") +
  scale_x_continuous(breaks = 0:12) +
  theme_minimal() +
  theme(legend.position = "bottom")



ggsave(
  filename="figures/diabetes_casestudy_trajectories.pdf",
  plot = trajectories,
  scale=1,
  width=10,
  height=6,
  units="in",
  dpi=600
)

## Table with analysis output for each method

library(purrr)

analysis_methods <- list(
  ## Treatment policy estimands
  "ipwtp" = analyse_diabetes_ipw(strategy = "treatment_policy"),
  "mmrmtp" = analyse_diabetes_mmrm(strategy = "treatment_policy"),
  "mitp" = analyse_diabetes_mi(strategy = "treatment_policy"),
  ## Hypothetical estimands
  "ipwhyp" = analyse_diabetes_ipw(strategy = "hypothetical"),
  "dmhyp" = analyse_diabetes_demediation(),
  "gcomhyp" = analyse_diabetes_gcomputation(),
  "mmrmhyp" = analyse_diabetes_mmrm(strategy = "hypothetical"),
  "mihyp" = analyse_diabetes_mi(strategy = "hypothetical")
)

set.seed(123)
analysis_tbl <-
  imap_dfr(
    analysis_methods,
    \(fit_fun, method_name) {
      fit <- fit_fun(cond, dat)
    }
  )

analysis_table <- analysis_tbl |> mutate(method = names(analysis_methods), .before = 1) |> select(c(method, "Estimate" = coef,
                                                                                                    "Standard error" = se,
                                                                                                    "95% CI lower" = ci_lower,
                                                                                                    "95% CI upper" = ci_upper,
                                                                                                    "p-value" = p))

knitr::kable(
  analysis_table,
  format = "latex",
  caption = "Table with analysis output for each method",
  label = "diabetes_casestudy_analysis"
) |>
  cat(file = "tables/diabetes_casestudy_analysis.tex")

