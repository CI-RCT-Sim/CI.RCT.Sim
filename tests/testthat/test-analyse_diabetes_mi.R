test_that("mi diabetes works", {
  Design <- assumptions_diabetes_rescue(print = FALSE) |>
    true_summary_statistics_diabetes_rescue()

  my_analyse <- analyse_diabetes_rescue_mi(
    ci_level = 0.95,
    estimand = "treatment_policy"
  )

  # Generate data with no missingness to check equivalence with linear regression
  withr::with_seed(123, {
    dat <- generate_diabetes_rescue(Design[7, ])
  })
  expect_no_error({
    res <- my_analyse(Design[7, ], dat)
  })

  mod_lm <- lm(y12 ~ trt + age + y0, data = dat)
  expect_equal(res$coef, coef(mod_lm)["trt"], tolerance = 1e-8)

  # Generate data with no missingness, and no treatment effect, to check equivalence with linear regression
  withr::with_seed(129, {
    dat <- generate_diabetes_rescue(Design[12, ])
  })
  expect_no_error({
    res <- my_analyse(Design[12, ], dat)
  })

  mod_lm <- lm(y12 ~ trt + age + y0, data = dat)
  expect_equal(res$coef, coef(mod_lm)["trt"], tolerance = 1e-8)

  # sanity checks under the null hypothesis
  expect_gt(res$p, 0.025)
  expect_lt(res$ci_lower, 0)
  expect_gt(res$ci_upper, 0)


  # Hypothetical strategy
  my_analyse <- analyse_diabetes_rescue_mi(
    ci_level = 0.95,
    estimand = "hypothetical"
  )

  # Generate data with no missingness to check equivalence with linear regression
  withr::with_seed(123, {
    dat <- generate_diabetes_rescue(Design[7, ])
  })
  expect_no_error({
    suppressWarnings(
      res <- my_analyse(Design[7, ], dat)
    )
  })
})
