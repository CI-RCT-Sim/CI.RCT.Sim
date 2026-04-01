#' Create Analyse Functions for
#'
#' @param recensor logical, whether to recensor the data after adjusting for treatment switching. Default is TRUE.
#' @param alpha numeric, significance level for confidence intervals. Default is 0.05.
#' @param B integer, number of bootstrap samples for estimating standard errors. Default is 100.
#'
#' @return an analyse function that can be used in runSimulation
#' @export
#'
#' @importFrom trtswitch tsesimp
#'
#' @examples
#' setting <- oncology_scenario()[1, ]
#' dat <- generate_oncology(setting)
#' analyse_oncology_TSE()(setting, dat)
analyse_oncology_TSE <- function(recensor = TRUE, alpha = 0.05, B = 100) {
  function(condition, dat, fixed_objects = NULL) {
    prep_data_RPSFTM_fun <- function(data) {
      data$time_on_trt <- 0
      data$time_on_trt[data$trt == 1] <- data$event_time[data$trt == 1]
      set_sw <- data$trt == 0 & data$switch == 1
      data$time_on_trt[set_sw] <- data$event_time[set_sw] - data$prog_time[set_sw]
      data$time_on_trt_relative <- data$time_on_trt / data$event_time
      data$max_FU <- data$calendar_end_of_study - data$calendar_start_time
      # for admin censoring, max_FU is the event time. The numeric calclation may cause minimal differences
      # so
      set_admin_cens <- data$max_FU < data$event_time
      data$max_FU[set_admin_cens] <- data$event_time[set_admin_cens]
      data
    }

    data_rpsftm <- prep_data_RPSFTM_fun(dat)

    TSE <- tsesimp(
      data = data_rpsftm,
      id = "id",
      time = "event_time",
      event = "ev",
      treat = "trt",
      censor_time = "max_FU",
      pd = "prog_ev",
      pd_time = "prog_time",
      swtrt = "switch",
      swtrt_time = "prog_time",
      base_cov = c("X_0", "W_0"),
      base2_cov = c("X_2BL", "W_2BL"),
      aft_dist = "weibull",
      alpha = alpha,
      recensor = recensor,
      swtrt_control_only = TRUE,
      offset = 0,
      boot = B > 0,
      n_boot = B
    )

    SE <- stats::sd(log(TSE$hr_boots))
    list(
      HR = TSE$hr,
      SElogHR = SE,
      low = TSE$hr_CI[1],
      up = TSE$hr_CI[2],
      p = TSE$pvalue,
      N_pat = nrow(dat),
      N_evt = sum(dat$ev)
    )
  }
}
