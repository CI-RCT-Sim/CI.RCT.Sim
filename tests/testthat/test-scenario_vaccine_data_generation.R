test_that("data generation vaccine scenario", {
  Design_test <- params_scenarios_grid(
    p_V = c(0.1, 0.3, 0), # probability for binary covariate prognostic for ICE and infection risk
    p_W = c(0.1, 0.3, 0), # probability for binary covariate prognostic for ICE and infection risk and modifying treatment effect
    lambda_post = -log(1-1/10)/365.25, # force of infection (baseline infection hazard) after 14 days, yearly incidence of 1/10
    gamma_0  = c(1), # choosen arbitrarily for testing
    gamma_W  = c(-0.8, 0), # regression parameters for compliance
    gamma_V  = c(0.5, 0),
    gamma_A  = c(-0.357, 0),
    gamma_AW = c(-0.3, 0),
    beta_V  = c(log(1.5), 0), # regression parameters for time to infection
    beta_W  = c(log(1.2), 0),
    beta_A1  = log(1-c(0, 0.5, 0.8, 0.9)), # choosen according to beta_A2
    beta_A2  = log(1-c(0, 0.8)), #choosen arbitrarily for testing
    beta_AW = c(log(0.8), 0),
    n_trt     = c(50), # study design parameters
    n_ctrl    = c(50),
    follow_up = c(365)
  )

  expect_no_error(generate_vaccine(Design_test[1,]))

})
