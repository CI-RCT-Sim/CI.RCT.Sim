test_that("analyse vaccine per protocol works", {
  Design <- vaccine_scenario() |>
    vaccine_scenario_set_gamma_0() |>
    vaccine_scenario_set_true_eff() |>
    vaccine_scenario_set_samplesize()

  my_analyse <- analyse_vaccine_pp(ci_level=0.95)

  # pV > 0
  # set seed to avoid the off-chance that none of the V are 1
  withr::with_seed(123, {
    dat <- generate_vaccine(Design[3, ])
  })
  expect_no_error({res <- my_analyse(Design[3, ], dat)})

  # pW > 0
  # set seed to avoid the off chance that none of the W are 1
  withr::with_seed(456, {
    dat2 <- generate_vaccine(Design[4, ])
  })

  expect_no_error({res2 <- my_analyse(Design[4, ], dat2)})

  # no effect
  # set seed to ensure that we are in a situation in which the test does not commit a type I error
  withr::with_seed(789, {
    condition3 <- Design |>
      subset(beta_A2==0) |>
      head(1)

    dat3 <- condition3 |>
      generate_vaccine()
  })

  # some sanity checks on p values and CIs
  expect_no_error({res3 <- my_analyse(condition3, dat3)})
  expect_gt(res3$p, 0.025)
  expect_lt(res3$VE_lower, 0.3)

  # no covariates
  withr::with_seed(123, {
    condition4 <- Design |>
      subset(p_V==0) |>
      subset(p_W==0) |>
      head(1)

    dat4 <- condition4 |>
      generate_vaccine()
  })

  expect_no_error({res4 <- my_analyse(condition4, dat4)})
})

test_that("per protocol analysis does not use unobserved covariates", {
  Design <- vaccine_scenario(print=FALSE) |>
    vaccine_scenario_set_gamma_0() |>
    vaccine_scenario_set_true_eff() |>
    vaccine_scenario_set_samplesize()

  my_analyse <- analyse_vaccine_pp(ci_level=0.95)

  condition <- Design[sample(1:nrow(Design), 1), ]

  withr::with_seed(123, {
    dat1 <- generate_vaccine(condition, fixed_objects = list(include_unobserved=FALSE))
  })

  withr::with_seed(123, {
    dat2 <- generate_vaccine(condition, fixed_objects = list(include_unobserved=TRUE))
  })

  res1 <- my_analyse(condition, dat1)
  res2 <- my_analyse(condition, dat2)

  expect_identical(res1, res2)
})
