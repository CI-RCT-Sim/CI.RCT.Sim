test_that("setting sample size works", {
  Design1 <- vaccine_scenario() |>
    vaccine_scenario_set_gamma_0() |>
    vaccine_scenario_set_true_eff() |>
    vaccine_scenario_set_samplesize()

  # finite sample sizes
  expect_true(all(is.finite(Design1$n_ctrl)))
  expect_true(all(is.finite(Design1$n_trt)))

  Design2 <- vaccine_scenario() |>
    vaccine_scenario_set_gamma_0() |>
    vaccine_scenario_set_true_eff() |>
    vaccine_scenario_set_samplesize(VE_H1 = 0.6)

  # larger sample sizes with smaller effect
  expect_true(all(Design2$n_ctrl > Design1$n_ctrl))

  Design3 <- vaccine_scenario() |>
    vaccine_scenario_set_gamma_0() |>
    vaccine_scenario_set_true_eff() |>
    vaccine_scenario_set_samplesize(alpha = 0.05)

  # smaller sample sizes with larger alpha
  expect_true(all(Design3$n_ctrl < Design1$n_ctrl))

  Design4 <- vaccine_scenario() |>
    vaccine_scenario_set_gamma_0() |>
    vaccine_scenario_set_true_eff() |>
    vaccine_scenario_set_samplesize(CSE=0.4)

  # larger sample sizes with larger clinically significant effect
  expect_true(all(Design4$n_ctrl > Design1$n_ctrl))

  Design5 <- vaccine_scenario() |>
    vaccine_scenario_set_gamma_0() |>
    vaccine_scenario_set_true_eff() |>
    vaccine_scenario_set_samplesize(power=0.9)

  # larger sample sizes with larger power
  expect_true(all(Design5$n_ctrl > Design1$n_ctrl))
})


test_that("two functions to generate parameters output the same parameters", {
  x <- vaccine_scenario()
  y <- vaccine_scenario_grid()

  expect_identical(x, y)
})
