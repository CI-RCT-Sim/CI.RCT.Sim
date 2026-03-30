## Set up a design
## Considering we have
sim_parameters <- vaccine_scenario_tweak(p_V = 0.5,
                                         overall_compliance = 0.5,
                                         gamma_V = 0.5,
                                         beta_V = log(c(1,1.5,0.5)),
                                         beta_A2 = 0) |>
  # mutate(overall_compliance = 0.5,
  #        dose_interval = 14,
  #        lambda_post = -log(1-(1/300))/(30/7)) |>
  vaccine_scenario_set_gamma_0() |>
  vaccine_scenario_set_true_eff() |>
  vaccine_scenario_set_samplesize() |>
  within({
    VE = 1-rr_ps
  })

## Check true VE
sim_parameters |> mutate(across(starts_with('beta'),~1-exp(.)))
## this is alright



sim_parameters |> mutate(across(starts_with('beta'),~1-exp(.))) -> tpar


## Simulate a large dataset (1e6 patients) runs some time!!
df <- purrr::map(1:nrow(sim_parameters),\(i){
  sim_parameters[i,] |>
    transform(n_trt = 1e6,
              n_ctrl=1e6) |>
    generate_vaccine()
})


## merge and join results with parameters
results <- df |> bind_rows(.id = 'scenario') |> as_tibble()
results <- tpar |> mutate(scenario = as.character(1:3)) |> dplyr::left_join(results,by='scenario')

## true cumulative incidence according to parameter setup
(1-exp(-(exp(cbind(1,sim_parameters$beta_V)) %*% diag(c(0, tpar$lambda_post[1])))*24))

## observed cumulative incidence (if you use !=0 VE group also by trt)
results |> group_by(scenario,V) |>
  dplyr::summarise(icd = mean(evt))
## cumulative incidence stimmt



