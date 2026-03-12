test_that("mmrm diabetes works", {
  Design <- assumptions_diabetes_rescue(print = FALSE) |>
    true_summary_statistics_diabetes_rescue()

  my_analyse <- analyse_diabetes_rescue_mmrm(
    ci_level = 0.95,
    strategy = "treatment_policy"
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
    dat <- generate_diabetes_rescue(Design[15, ])
  })
  expect_no_error({
    res <- my_analyse(Design[15, ], dat)
  })

  mod_lm <- lm(y12 ~ trt + age + y0, data = dat)
  expect_equal(res$coef, coef(mod_lm)["trt"], tolerance = 1e-8)

  # sanity checks under the null hypothesis
  expect_gt(res$p, 0.025)
  expect_lt(res$ci_lower, 0)
  expect_gt(res$ci_upper, 0)
})
