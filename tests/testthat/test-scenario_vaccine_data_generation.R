test_that("data generation vaccine scenario", {
  Design_test <- vaccine_scenario() |>
    vaccine_scenario_set_gamma_0()

  Design_test$beta_A1 <- log(1-0.8)

  expect_no_error(generate_vaccine(Design_test[1,]))

})
