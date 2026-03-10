test_that("multiplication works", {
  Design <- vaccine_scenario() |>
    vaccine_scenario_set_gamma_0() |>
    vaccine_scenario_set_true_eff() |>
    vaccine_scenario_set_samplesize()

  condition <- Design[3, ]

  # combinations of non-estimable variables
  cond1 <- cond2 <- cond3 <- condition
  cond2$p_W <- 0.3
  cond2$p_V <- 0
  cond3$p_W <- 0.3

  dat1 <- generate_vaccine(cond1)
  dat2 <- generate_vaccine(cond2)
  dat3 <- generate_vaccine(cond3)

  my_analyse <- analyse_vaccine_pp(ci_level=0.95, VE_margin=0.3)

  expect_no_error({res1 <- my_analyse(cond1, dat1)})
  expect_no_error({res2 <- my_analyse(cond2, dat2)})
  expect_no_error({res3 <- my_analyse(cond3, dat3)})

  res1
  res1b <- analyse_vaccine_pp(VE_margin=0.2)(cond1, dat1)
  expect_lte(res1b$p, res1$p, "expect lower p-value for lower super-superiority margin")
})
