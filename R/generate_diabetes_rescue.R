#' Generate Dataset ...
#'
#' @param condition condition row of Design dataset
#' @param fixed_objects fixed objects of Design dataset
#'
#' @details
#' Condition has to contain the following columns:
#'
#'   * eff the effect size
#'   *
#'   * ...
#'
#' @return
#' For generate_diabetes_rescue: A data set with the columns id, trt
#' (1=treatment, 0=control), evt (event, currently TRUE for all observations)
#'
#' @export
#' @describeIn generate_diabetes_rescue simulates a dataset with ...
#'
#' @examples
generate_diabetes_rescue <- function(condition, fixed_objects=NULL){

  # What still remains to be done is to
    # incorporate sample size calculations
    # implement the correlation rho between the visits (which also goes in to the sample size calculations)
    # time varying treatment effect (and rescue effect)
    # documentation
    # Statements to be able to turn rescue and missingness on and off

  if (length(unique(condition$rescue_effect))!=condition$k[1]+1){
    stop(gettext("Non-matching dimensions"))
  } else if (length(unique(condition$eff))!=condition$k[1]+1){
    stop(gettext("Non-matching dimensions"))
  }

  n <- 100 # Should be edited according to proposed power and effect size
  # n <- ((qnorm(1 - alpha / 2) + qnorm(power))^2)*sigma^2/(eff^2)

  visit <- 0:condition$k[1]
  id=1:n
  trt=rbinom(n, 1, 0.5)
  age = rnorm(n, mean = condition$mean_age, sd = condition$sd_age)
  age_slope = 2 * exp(-condition$b_age[1] * (age - 30))
  response_trt = runif(n)
  mu <- matrix(NA, nrow = n, ncol = length(visit))
  for (i in 1:length(visit)){
    mu[,i] <- condition$mean_bl[1] +
      visit[i] / condition$k[1] * age_slope +
      unique(condition$eff)[i] * response_trt * trt
  }

  resid <- mvtnorm::rmvnorm(n, rep(0, length(visit)), diag(length(visit)))
  Y <- mu + resid

  h_0 <- log(condition$pr_rescue[1] / (1 - condition$pr_rescue[1]))
  expit <- function(x) exp(x) / (1 + exp(x))
  pr_rescue <- expit(h_0 + (Y - 10) * condition$h_y[1] + (age - condition$mean_age[1]) * condition$h_age[1])
  pr_rescue[,1] <- 0
  pr_rescue[,condition$k+1] <- 0

  resc <- matrix(rbinom((condition$k[1] + 1)*n, size = 1, prob = pr_rescue), nrow = n)
  rescue <- t(apply(resc, 1, cumsum))>0
  rescue_start <- rowSums(!rescue)+1
  k_rescue <- rowSums(rescue)
  response_rescue <- runif(n)
  any_rescue <- c()

  for (i in 1:n){
    if (k_rescue[i] > 0){
      rescue_set <- (rescue_start[i] + 1):(condition$k[1] + 1)
      # browser()
      Y[i,rescue_set] <- mu[i,rescue_set] +
        response_rescue[i] * unique(condition$rescue_effect)[rescue_set - rescue_start[i] + 1] +
        resid[rescue_set]
      any_rescue[i] <- TRUE
    } else{
      rescue_start[i] <- NA
      any_rescue[i] <- FALSE
    }
  }

  g_0 <- log(condition$pr_missing[1] / (1 - condition$pr_missing[1]))
  pr_miss <- expit(g_0 + (Y - 10) * condition$g_y[1] +
                     (age - condition$mean_age[1]) * condition$g_age[1] +
                     rescue * condition$g_rescue[1]) # actual prob. to drop out
  pr_miss[,1] <- 0 # we assume complete data at baseline

  wd <- matrix(rbinom((condition$k[1] + 1)*n, size = 1, prob = pr_miss), nrow = n)
  wd1 <- t(apply(wd, 1, cumsum))>0

  for (i in 1:n){
    miss_start <- sum(!wd1[i,]) + 1
    if (miss_start <= (condition$k[1] + 1)) Y[i,miss_start:(condition$k[1] + 1)] <- NA
  }

  out <- data.frame(id, trt, age, Y, any_rescue, rescue_start)
  names(out) <- c("id", "trt", "age", paste("y", visit, sep = ""), "any_rescue", "rescue_start")
  out
}

#' Create an empty assumptions data.frame for generate_diabetes_rescue
#'
#' @param print print code to generate parameter set?
#'
#' @return For assumptions_diabetes_rescue: a design tibble with default values invisibly
#'
#' @details assumptions_diabetes_rescue generates a default design `data.frame`
#'   for use with generate_diabetes_rescue If print is `TRUE` code to produce
#'   the template is also printed for copying, pasting and editing by the user.
#'   (This is the default when run in an interactive session.)
#'
#' @export
#' @describeIn generate_diabetes_rescue generate default design tibble
#'
#' @examples
#' Design <- assumptions_diabetes_rescue()
#' Design
assumptions_diabetes_rescue <- function(print=interactive()){
  skel <- "expand.grid(
  eff = c(0,1,2,3,4,5),                  # treatment effect
  rescue_effect = c(0,-2,-4,-6,-8,-10),  # effect of rescue medication
  k = 5,                                 # Number of visits post baseline
  mean_bl=8,                             # mean hbalc value at baseline
  mean_age=60,                           # mean of the variable age
  sd_age=10,                             # standard deviation of the variable age
  b_age=log(2)/10,                       # age coefficient
  pr_rescue = 0.05,                      # probability for rescue medication
  h_y = log(3),                          # strong effect due to high hba1c
  h_age = -log(1.01),                    # weaker age effect than for dropout
  pr_missing = 0.02,                     # probability for missing data
  g_y = log(1.5),                        # moderate effect due to high hba1c
  g_age = log(1.02),                     # older patients drop out more easily, let's say, stronger age effect than for rescue
  g_rescue = log(1.5)                    # notable effect due to rescue medication, increase to large value to have positivity violation, like 5 or 10
)
"

if(print){
  cat(skel)
}

invisible(
  skel |>
    str2expression() |>
    eval()
)
}

#' Calculate true summary statistics for scenarios with delayed treatment effect
#'
#' @param Design Design data.frame for x
#' @param cutoff_stats Cutoff time for rmst and average hazard ratios
#' @param fixed_objects fixed objects not used for now
#'
#' @return For true_summary_statistics_x: the design data.frame
#'   passed as argument with the additional columns:
#' * `rmst_trt` rmst in the treatment group
#' * `median_surv_trt` median survival in the treatment group
#' * `rmst_ctrl` rmst in the control group
#' * `median_surv_ctrl` median survial in the control group
#' * `gAHR` geometric average hazard ratio
#' * `AHR` average hazard ratio
#'
#' @export
#'
#' @describeIn generate_diabetes_rescue  calculate true summary statistics for ...
#'
#' @examples
true_summary_statistics_diabetes_rescue <- function(Design, cutoff_stats=10, fixed_objects=NULL){

  true_summary_statistics_diabetes_rescue_rowwise <- function(condition, cutoff_stats){
    res <- data.frame(
      rmst_trt = NA_real_,
      medial_surv_trt = NA_real_,
      rmst_ctrl = NA_real_,
      median_surv_ctrl = NA_real_,
      gAHR = NA_real_,
      AHR = NA_real_
    )

    res
  }

  Design <- Design |>
    split(1:nrow(Design)) |>
    mapply(FUN=true_summary_statistics_x_rowwise, cutoff_stats = cutoff_stats, SIMPLIFY = FALSE)

  Design <- do.call(rbind, Design)

  Design
}
