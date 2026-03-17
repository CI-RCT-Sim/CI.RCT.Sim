test_that("demediation diabetes works", {
  Design <- diabetes_scenario(print = FALSE) |>
    diabetes_scenario_set_truevalues()

  my_analyse <- analyse_diabetes_demediation()

  # Generate data
  withr::with_seed(123, {
    dat <- generate_diabetes(Design[1, ])
  })
  expect_no_error({
    res <- my_analyse(Design[1, ], dat)
  })

  # Generate data with no missingness, and no treatment effect
  withr::with_seed(129, {
    dat <- generate_diabetes(Design[15, ])
  })
  expect_no_error({
    res <- my_analyse(Design[15, ], dat)
  })

  # sanity checks under the null hypothesis
  expect_lt(res$ci[1], 0)
  expect_gt(res$ci[2], 0)
})
