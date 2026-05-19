#' Create Analyse Functions for ITT analysis
#'
#' @param X input can be used to pass parameters to the analyse function
#'
#' @return an analyse function that can be used in runSimulation
#' @export
#'
#' @importFrom survival coxph Surv
#' @importFrom stats confint
#'
#' @examples
#' setting <- oncology_scenario()[1, ]
#' data <- generate_oncology(setting)
#' analyse_oncology_itt()(setting, data)
analyse_oncology_itt <- function(X) {
  function(condition, data, fixed_objects = NULL) {
    mod <- coxph(Surv(time = event_time, event = ev) ~ trt + X_0 + W_0, data = data)
    HR <- exp(coef(mod)[1])
    CI <- exp(confint(mod)[1, ])
    smr <- summary(mod)
    p <- smr$coef[1, "Pr(>|z|)"]
    p<-ifelse(HR<=1,p/2,1-p/2)
    SE <- smr$coef[1, "se(coef)"]
    list(
      HR = HR,
      SElogHR = SE,
      low = CI[[1]],
      up = CI[[2]],
      p = p,
      N_pat = nrow(data),
      N_evt = sum(data$ev)
    )
  }
}
