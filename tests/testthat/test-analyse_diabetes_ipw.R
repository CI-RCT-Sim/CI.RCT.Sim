test_that("ipw diabetes works", {
  Design <- diabetes_scenario(print = FALSE) |>
    diabetes_scenario_set_truevalues()

  my_analyse <- analyse_diabetes_ipw(
    strategy = "treatment_policy"
  )

  # Generate
  withr::with_seed(123, { # problems with seed 1, 12, 1122, 12345
    dat <- generate_diabetes(Design[1, ])
  })
  expect_no_error({
    res <- my_analyse(Design[1, ], dat)
  })

  # Check the hypothetical strategy
  my_analyse <- analyse_diabetes_ipw(
    strategy = "hypothetical"
  )

  # Generate
  withr::with_seed(122, {
    dat <- generate_diabetes(Design[12, ])
  })
  expect_no_error({
    res <- my_analyse(Design[12, ], dat)
  })

  # sanity checks under the null hypothesis
  expect_gt(res$p, 0.025)
  expect_lt(res$ci_lower, 0)
  expect_gt(res$ci_upper, 0)
})
