test_that("generate oncology works", {
  Design <- oncology_scenario() |>
    oncology_scenario_set_truevalues()

  # Generate data
  expect_no_error({
    withr::with_seed(123, {
      dat <- generate_oncology(Design[1, ])
    })
  })
})
