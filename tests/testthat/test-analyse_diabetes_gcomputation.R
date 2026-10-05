test_that("g computation diabetes works", {
  Design <- diabetes_scenario(print = FALSE) |>
    diabetes_scenario_set_truevalues()

  my_analyse <- analyse_diabetes_gcomputation()

  # Generate data with no missingness to check equivalence with linear regression
  withr::with_seed(123, {
    dat <- generate_diabetes(Design[7, ])
  })
  expect_no_error({
    res <- my_analyse(Design[7, ], dat)
  })

  # Generate data with no missingness, and no treatment effect
  withr::with_seed(12, {
    dat <- generate_diabetes(Design[9, ])
  })
  expect_no_error({
    res <- my_analyse(Design[9, ], dat)
  })

  # sanity checks under the null hypothesis
  expect_lt(res$ci_lower, 0)
  expect_gt(res$ci_upper, 0)
})
