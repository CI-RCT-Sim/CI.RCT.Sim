
test_that("ivreg vaccine works", {
  Design <- vaccine_scenario() |>
    vaccine_scenario_set_gamma_0() |>
    vaccine_scenario_set_true_eff() |>
    vaccine_scenario_set_samplesize()

  my_analyse <- analyse_vaccine_ivreg(ci_level=0.95)

  # pV > 0

  dat <- generate_vaccine(Design[3, ])
  expect_no_error({res <- my_analyse(Design[3, ], dat)})

  dat1 <- dat |>
    within({
      T <- factor(ifelse(((trt==1) & (C==1)), 1L, 0L), levels=c("1", "0"))
      V <- factor(V)
      trt <- factor(trt)
      evt <- as.integer(evt)
    })
  mod_ivreg <- ivreg::ivreg(evt ~ T + V | trt + V, data=dat1)

  expect_lt(abs(mod_ivreg$coefficients["T0"] + res$RD), .Machine$double.eps * 2, "same result as with IV reg from package")

  # pW > 0

  dat2 <- generate_vaccine(Design[4, ])

  expect_no_error({res2 <- my_analyse(Design[4, ], dat2)})

  dat1 <- dat2 |>
    within({
      T <- factor(ifelse(((trt==1) & (C==1)), 1L, 0L), levels=c("1", "0"))
      W <- factor(W)
      trt <- factor(trt)
      evt <- as.integer(evt)
    })
  mod_ivreg2 <- ivreg::ivreg(evt ~ T + W | trt + W, data=dat1)

  expect_lt(abs(mod_ivreg2$coefficients["T0"] + res2$RD), .Machine$double.eps * 2, "same result as with IV reg from package")

})
