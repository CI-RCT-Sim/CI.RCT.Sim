devtools::load_all()

Design <- params_scenarios_grid(
  p_V = c(0.1, 0.3), # probability for binary covariate prognostic for ICE and infection risk
  p_W = c(0.1, 0.3), # probability for binary covariate prognostic for ICE and infection risk and modifying treatment effect
  lambda_post = 100* -log(1-(1/c(500, 1000, 2000)))/365.25, # force of infection (baseline infection hazard) after 14 days, yearly incidence of 1/500, 1/1000, 1/2000
  overall_compliance = c(0.95), # used to callibrate gamma0
  gamma_W  = c(-0.8, 0), # regression parameters for compliance
  gamma_V  = c(0.5, 0),
  gamma_A  = c(-0.357, 0),
  gamma_AW = c(-0.3, 0),
  beta_V  = c(log(1.5), 0), # regression parameters for time to infection
  beta_W  = c(log(1.2), 0),
  effect_before_d2 = c(1,0), # indicator whether there's any effect before d2, used to set beta_A1
  beta_A2  = log(1-c(0.8, 0, 0.5, 0.9)), #
  beta_AW = c(log(0.8), 0),
  n_trt     = c(500), # study design parameters
  n_ctrl    = c(500),
  follow_up = c(365)
) |>
  transform(
    # either 0 (no effect before d2) or VE before dose 2 is VE after dose 2 - 0.3
    beta_A1 = log(1-pmax(1-exp(beta_A2)-0.3, 0.))*effect_before_d2
  )

Design <- Design |>
  vaccine_scenario_set_gamma_0() |>
  vaccine_scenario_set_true_eff()

condition <- Design[7, ]
n_sim <- 10000

test <- replicate(n_sim, {
  sample <- generate_vaccine(condition, fixed_objects = list(include_unobserved=TRUE))

  ps <- sample |>
    subset(C_1==1) |>
    subset(C_0==1)

  risks <- ps |>
    by(~trt, \(grp){
      mean(grp$evt)
    })


    rr <- (risks[["1"]] / risks[["0"]])
    rr
})

hist(test)
abline(v=condition$rr_ps, lwd=1.5, col="blue")

# significant difference between calculation and simulation?
t.test(test, mu=condition$rr_ps)

# close, but still quite a bit off with high confidence
