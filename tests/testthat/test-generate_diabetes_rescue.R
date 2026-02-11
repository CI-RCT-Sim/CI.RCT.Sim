test_that("assumptions_diabetes_rescue works", {
  assumptions <- assumptions_diabetes_rescue()

  expect_named(assumptions, c(
    "k", "mean_age", "sd_age", "b_age", "mean_bl", "sd_bl", "rho", "delta",
    "lambda", "delta_resc", "lambda_resc", "resc_0", "resc_y", "resc_age",
    "miss_0", "miss_y", "miss_age", "miss_resc", "hyp"
  ))

  expect_equal(ncol(assumptions), 19)
  expect_equal(nrow(assumptions), 22)
  # whole data set should be duplicated for values 0 and 1 of hyp
  expect_equal(assumptions$hyp, c(rep(0, 11), rep(1, 11)))
  expect_equal(assumptions[1:11, -19], assumptions[12:22, -19], ignore_attr=TRUE)
  # rows should only differ in one variable for the first and the second half
  # (reference scenarion, then differ in one parameter for each row)
  expect_equal((assumptions[rep(1, 10),] != assumptions[2:11,]) |> rowSums(), rep(1, 10), ignore_attr=TRUE)
  expect_equal((assumptions[rep(12, 10),] != assumptions[13:22,]) |> rowSums(), rep(1, 10), ignore_attr=TRUE)
})



