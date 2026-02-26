test_that("true summary stats for minimal example work", {
  res <- true_summary_statistics_minimal_example(data.frame(mean1=2, mean0=1))
  expect_equal(res$eff_size, 1.)
})

test_that("assumptions minimal example work", {

  expected_output <- tibble::tribble(
    ~n, ~mean1, ~mean0,
    50,      1,      0,
   100,      1,      0,
    50,      0,      0,
  )
  assumptions <- expect_output(assumptions_minimal_example(print=TRUE))
  expect_equal(assumptions, expected_output)

})

test_that("generating minimal example works", {
  output <- generate_minimal_example(data.frame(n=10, mean1=0, mean0=0), fixed_objects = NULL)
  expect_named(output, c("group", "y"))
  expect_equal(nrow(output), 20)
  expect_type(output$y, "double")
  expect_true(all(output$group %in% c(0.,1.)))
})
