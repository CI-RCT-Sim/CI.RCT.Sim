#' Create Analyse Functions for IPW
#'
#' @param X input can be used to pass parameters to the analyse function
#'
#' @return an analyse function that can be used in runSimulation
#' @export
#'
#' @importFrom survival tmerge coxph Surv
#' @importFrom stats confint glm predict
#'
#' @examples
#' setting <- oncology_scenario()[1, ]
#' dat <- generate_oncology(setting)
#' analyse_oncology_ipw()(setting, dat)
analyse_oncology_ipw <- function(X) {
  function(condition, dat, fixed_objects = NULL) {
    set <- dat$trt == 0 & dat$prog_ev == 1
    wmod <- glm(switch ~ X_2BL + W_2BL, family = binomial, data = dat, subset = set)
    pred <- predict(wmod, type = "response")
    pred_a <- ifelse(dat$switch[set] == 1, pred, 1 - pred)
    dat$w <- 1
    dat$w[set] <- 1 / pred_a

    sdat <- tmerge(data1 = dat, data2 = dat, id = id, tstop = event_time)
    sdat <- tmerge(data1 = sdat, data2 = dat, id = id, death_event = event(event_time, ev))
    sdat <- tmerge(data1 = sdat, data2 = dat, id = id, PD = tdc(prog_time))
    sdat$w[sdat$PD == 0 | sdat$trt == 1] <- 1

    remo <- sdat$trt == 0 & sdat$PD == 1 & sdat$switch == 1
    dat_a <- sdat[!remo, ]
    mod <- coxph(
      Surv(time = tstart, time2 = tstop, event = death_event) ~ trt + X_0 + W_0,
      data = dat_a,
      weights = w,
      robust = TRUE,
      id = id
    )

    HR <- exp(coef(mod)[1])
    CI <- exp(confint(mod)[1, ])
    smr <- summary(mod)
    p <- smr$coef[1, "Pr(>|z|)"]
    SE <- smr$coef[1, "robust se"]
    list(
      HR = HR,
      SElogHR = SE,
      low = CI[1],
      up = CI[2],
      p = p,
      N_pat = nrow(dat),
      N_evt = sum(dat$ev)
    )
  }
}
