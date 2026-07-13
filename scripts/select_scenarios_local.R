library(dplyr)
library(glue)

source("scripts/vaccine_scenario_classes.R")

all_results <- bind_rows(
  vaccine_scenario_A1 |> mutate(scen_cat="A1", row=1:n()),
  vaccine_scenario_A2 |> mutate(scen_cat="A2", row=1:n()),
  vaccine_scenario_B1 |> mutate(scen_cat="B1", row=1:n()),
  vaccine_scenario_C1 |> mutate(scen_cat="C1", row=1:n()),
  vaccine_scenario_D1 |> mutate(scen_cat="D1", row=1:n()),
  vaccine_scenario_extra |> mutate(scen_cat="extra", row=1:n()),
)|>
  vaccine_scenario_set_beta_A1_relative() |>
  vaccine_scenario_set_gamma_0() |>
  vaccine_scenario_set_true_eff() |>
  vaccine_scenario_set_samplesize() |>
  within({
    VE = 1-rr_ps
    scenario_nr = seq_along(VE)
  })


all_results |>
  group_by(scen_cat) |>
  filter(n_trt == 7639,
         round(exp(beta_A2),1) %in% c(0.3,0.7)) |>
  filter((scen_cat == 'A1') |
           (scen_cat == 'A2' & gamma_A == -0.357) |
           (scen_cat == 'B1' & gamma_A == -0.357 & gamma_V == 0.5 & round(beta_V,1) == 0.4) |
           (scen_cat == 'C1' & gamma_A == -0.357 & gamma_W == -0.8 & gamma_AW == -0.3 & round(beta_W,2) == 0.18 & round(beta_AW,2) == -0.22) |
           (scen_cat == 'D1' & gamma_A == -0.357 & gamma_V == 0.5 & gamma_W == -0.8 & gamma_AW == -0.3 & round(beta_V,1) == 0.4 & round(beta_W,2) == 0.18 & round(beta_AW,2) == -0.22)
  ) |>
  select(scen_cat,row) |>
  rowwise() |>
  mutate(
    bash=glue("scenario={scen_cat}    SLURM_ARRAY_TASK_ID={row} RENV_CONFIG_SANDBOX_ENABLED=FALSE Rscript ./scripts/vaccine_rowwise.R")
  ) |>
  ungroup() |>
  pull(bash) |>
  paste(collapse="\n") |>
  cat()
