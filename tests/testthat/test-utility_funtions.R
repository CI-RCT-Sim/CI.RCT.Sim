test_that("generating grid of simulation parameters works", {
  example_value <- params_scenarios_grid(x=1:3, y=letters[1:3], z=factor("A"))
  expected_value <- tibble::tibble(
    x=c(1,2,3,1,1),
    y=letters[c(1,1,1,2,3)],
    z=factor("A")
  )

  expect_equal(example_value, expected_value)

  example_value_list <- params_scenarios_grid(x=1:2, y=list(c(1), c(1, 2)))
  expected_value_list <- tibble::tibble(
    x=c(1,2,1),
    y=list(c(1), c(1), c(1,2))
  )

  expect_equal(example_value_list, expected_value_list)
})
