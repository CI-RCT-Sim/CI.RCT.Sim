#' Create Analyse Functions for
#'
#' @param X input can be used to pass parameters to the analyse function
#'
#' @return an analyse function that can be used in runSimulation
#' @export
#'
#' @examples
#' \donttest{
#' setting <- oncology_scenario()[1, ]
#' dat <- generate_oncology(setting)
#' analyse_oncology_cens()(setting, dat)
#' }
analyse_oncology_cens <- function(X) {
  function(condition, dat, fixed_objects = NULL) {
    set_cens <- data$trt == 0 & data$switch == 1
    data$event_time[set_cens] <- data$prog_time[set_cens]
    data$ev[set_cens] <- 0
    mod <- coxph(Surv(time = event_time, event = ev) ~ trt + X_0 + W_0, data = data)
    HR <- exp(coef(mod)[1])
    CI <- exp(confint(mod)[1, ])
    smr <- summary(mod)
    p <- smr$coef[1, "Pr(>|z|)"]
    SE <- smr$coef[1, "se(coef)"]
    list(
      HR = HR,
      SElogHR = SE,
      low = CI[1],
      up = CI[2],
      p = p
    )
  }
}
