test_that("generating grid of simulation parameters works", {
  example_value <- params_scenarios_grid(x=1:3, y=letters[1:3], z=factor("A"))
  expected_value <- data.frame(
    x=c(1,2,3,1,1),
    y=letters[c(1,1,1,2,3)],
    z=factor("A")
  )

  expect_equal(example_value, expected_value)
})
