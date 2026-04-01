test_that("TSE oncology works", {
  Design <- oncology_scenario(print = FALSE) |>
    oncology_scenario_set_truevalues()

  my_analyse <- analyse_oncology_TSE(recensor = TRUE)

  # Generate data
  withr::with_seed(123, {
    dat <- generate_oncology(Design[1, ])
  })
  expect_no_error({
    res <- my_analyse(Design[1, ], dat)
  })

  # Generate data with no missingness, and no treatment effect, to check equivalence with linear regression
  withr::with_seed(129, {
    dat <- generate_oncology(Design[50, ])
  })
  expect_no_error({
    res <- my_analyse(Design[50, ], dat)
  })

  # sanity checks under the null hypothesis
  expect_gt(res$p, 0.025)
  expect_lt(res$low, 1)
  expect_gt(res$up, 1)

  my_analyse <- analyse_oncology_TSE(recensor = FALSE)

  # Generate data
  withr::with_seed(123, {
    dat <- generate_oncology(Design[1, ])
  })
  expect_no_error({
    res <- my_analyse(Design[1, ], dat)
  })

  # Generate data with no missingness, and no treatment effect, to check equivalence with linear regression
  withr::with_seed(129, {
    dat <- generate_oncology(Design[50, ])
  })
  expect_no_error({
    res <- my_analyse(Design[50, ], dat)
  })

  # sanity checks under the null hypothesis
  expect_gt(res$p, 0.025)
  expect_lt(res$low, 1)
  expect_gt(res$up, 1)
})
