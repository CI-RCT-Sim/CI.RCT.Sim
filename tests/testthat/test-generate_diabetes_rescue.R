test_that("assumptions_diabetes_rescue works", {
  assumptions <- assumptions_diabetes_rescue()

  expect_named(assumptions, c(
    "k", "mean_age", "sd_age", "b_age", "mean_bl", "sd_bl", "rho", "delta",
    "lambda", "delta_resc", "lambda_resc", "resc_0", "resc_y", "resc_age",
    "miss", "hyp"
  ))

  expect_equal(ncol(assumptions), 16)
  expect_equal(nrow(assumptions), 14)
  # whole data set should be duplicated for values 0 and 1 of hyp
  expect_equal(assumptions$hyp, c(rep(1, 7), rep(0, 7)))
  expect_equal(assumptions[0:7, -16], assumptions[8:14, -16], ignore_attr = TRUE)
  # rows should only differ in one variable for the first and the second half
  # (reference scenarion, then differ in one parameter for each row)
  expect_equal((assumptions[rep(1, 4), -15] != assumptions[2:5, -15]) |> rowSums(), rep(1, 4), ignore_attr = TRUE)
  expect_equal((assumptions[rep(8, 4), -15] != assumptions[9:12, -15]) |> rowSums(), rep(1, 4), ignore_attr = TRUE)
})
