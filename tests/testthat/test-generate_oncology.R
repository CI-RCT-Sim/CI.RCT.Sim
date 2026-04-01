test_that("generate oncology works", {
  Design <- oncology_scenario(print = FALSE) |>
    oncology_scenario_set_truevalues()

  # Generate data
  withr::with_seed(123, {
    dat <- generate_oncology(Design[1, ])
  })
})
