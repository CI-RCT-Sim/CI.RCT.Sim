setwd("C:\\EMA_Causal\\CI.RCT.Sim")
source("scripts\\plots_common.R")

results <- read_multiple(
  dir_ls("Z:\\EMA_causal_16May_node2\\CI.RCT.Sim\\results2", regexp = "May_1K_scen_small_n_H1_*"),
  check_design_names = FALSE
)
results$rpsftm.test.rejection_0.025


# install.packages("readxl")
table <- readxl::read_xlsx("results/oncology_scenario_list.xlsx")
tab <- table# |> dplyr::filter(!grepl("Extra", Description))

# Order them correctly
results_long <- results |>
  mutate(
    onco_start = as.numeric(str_extract(file, "(?<=onco)\\d+"))
  ) |>
  arrange(onco_start) |>
  #mutate(
  #  scenario_nr = 1:n(),
  #  class = case_when(
  #    scenario_nr <= 24 ~ "higheff",
  #    scenario_nr <= 47 ~ "loweff",
  #    scenario_nr <= 69 ~ "higheff_null",
  #    scenario_nr <= 91 ~ "loweff_null"
  #  ),
  #  scenario_names = tab$scenario_name
  #) |>
  results_pivot_longer()


lapply(c("higheff", "loweff", "higheff_null", "loweff_null"), \(x){
  results_long <- results_long |>
    filter(class == x) |>
    mutate(scenario_names = factor(scenario_names, levels = unique(scenario_names)))

  pos <- position_jitterdodge(dodge.width = 0.8, jitter.width = 0.1)

  # Plots for Estimators ----------------------------------------------------

  plot_bias <- results_long |>
    filter(str_detect(method, "\\.est")) |>
    ggplot(aes(x = scenario_names, colour = method, group = method)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_colour_brewer(palette = "Set1") +
    aes(y = bias) +
    # geom_line() +
    geom_point(position = pos) +
    geom_errorbar(aes(ymin = bias - sd_bias / sqrt(REPLICATIONS), ymax = bias + sd_bias / sqrt(REPLICATIONS)), width = 0.5, position = pos) +
    geom_hline(yintercept = 0)

  save_plot(plot_bias, paste0("figures/onco_", x, "_bias.pdf"))

  plot_bias_restricted <- results_long |>
    filter(str_detect(method, "\\.est")) |>
    ggplot(aes(x = scenario_names, colour = method, group = method)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_colour_brewer(palette = "Set1") +
    aes(y = bias) +
    # geom_line() +
    geom_point(position = pos) +
    geom_errorbar(aes(ymin = bias - sd_bias / sqrt(REPLICATIONS), ymax = bias + sd_bias / sqrt(REPLICATIONS)), width = 0.5, position = pos) +
    geom_hline(yintercept = 0) +
    scale_y_continuous(limits = c(-0.5, 0.75))

  save_plot(plot_bias_restricted, paste0("figures/onco_", x, "_bias_restricted.pdf"))

  plot_coverage <- results_long |>
    filter(str_detect(method, "\\.est")) |>
    ggplot(aes(x = scenario_names, colour = method, group = method)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_colour_brewer(palette = "Set1") +
    aes(y = coverage) +
    # geom_line() +
    geom_point(position = pos) +
    geom_hline(yintercept = 0.95) +
    scale_y_continuous(limits = c(0, 1))

  save_plot(plot_coverage, paste0("figures/onco_", x, "_coverage.pdf"))

  plot_coverage_restricted <- results_long |>
    filter(str_detect(method, "\\.est")) |>
    ggplot(aes(x = scenario_names, colour = method, group = method)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_colour_brewer(palette = "Set1") +
    aes(y = coverage) +
    # geom_line() +
    geom_point(position = pos) +
    geom_hline(yintercept = 0.95) +
    scale_y_continuous(limits = c(0.5, 1))

  save_plot(plot_coverage_restricted, paste0("figures/onco_", x, "_coverage_restricted.pdf"))

  plot_ci_width <- results_long |>
    filter(str_detect(method, "\\.est")) |>
    ggplot(aes(x = scenario_names, colour = method, group = method)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_colour_brewer(palette = "Set1") +
    aes(y = width) +
    # geom_line() +
    geom_point(position = pos) +
    coord_cartesian(ylim = c(0, 4))

  save_plot(plot_ci_width, paste0("figures/onco_", x, "_ci_width.pdf"))

  # Plots for Tests ---------------------------------------------------------

  plot_rejection_tests <- results_long |>
    filter(str_detect(method, "\\.test")) |>
    ggplot(aes(x = scenario_names, colour = method, group = method)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_colour_brewer(palette = "Set1") +
    aes(y = rejection_0.05) +
    # geom_line() +
    geom_point(position = pos) +
    geom_hline(yintercept = 0.05) +
    scale_y_continuous(limits = c(0, 1))

  save_plot(plot_rejection_tests, paste0("figures/onco_", x, "_rejection_tests.pdf"))

  plot_rejection_est <- results_long |>
    filter(str_detect(method, "\\.est")) |>
    ggplot(aes(x = scenario_names, colour = method, group = method)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_colour_brewer(palette = "Set1") +
    aes(y = 1 - null_lower) +
    # geom_line() +
    geom_point(position = pos) +
    geom_hline(yintercept = 0.05) +
    scale_y_continuous(limits = c(0, 1))

  save_plot(plot_rejection_est, paste0("figures/onco_", x, "_rejection_est.pdf"))
})
