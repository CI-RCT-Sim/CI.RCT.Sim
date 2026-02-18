test_that("analyse minimal example works", {
  test_data <- data.frame(
    y = c(0.1, -0.1, 0.9, 1.1),
    group = c(0,0,1,1)
  )

  # results should be the same as for t-test
  t_test <- t.test(y~I(1-group), data=test_data)
  expected_lm <- list(
    p=t_test$p.value,
    coef=c(group=1),
    ci_lower=t_test$conf.int[1],
    ci_upper=t_test$conf.int[2]
  )

  result_lm <- analyse_minimal_example_lm()(condition=NA, dat=test_data)
  expect_equal(result_lm, expected_lm)

  expected_t <- list(
    p = 2*pt(1/sqrt(0.1^2 + 0.1^2), 2, lower.tail = FALSE)
  )
  result_t <- analyse_minimal_example_t()(condition=NA, dat=test_data)

  expect_equal(expected_t, result_t)
})
