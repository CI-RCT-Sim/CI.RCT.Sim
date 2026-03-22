#' Create Analyse Functions for
#'
#' @param X input can be used to pass parameters to the analyse function
#'
#' @return an analyse function that can be used in runSimulation
#' @export
#'
#' @importFrom trtswitch rpsftm
#'
#' @examples
#' \donttest{
#' setting <- oncology_scenario()[1, ]
#' dat <- generate_oncology(setting)
#' analyse_oncology_rpsftm()(setting, dat)
#' }
analyse_oncology_rpsftm <- function(recensor = TRUE, alpha = 0.05, B = 100) {
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

    RPS <- rpsftm(
      data = data_rpsftm, id = "id", time = "event_time", event = "ev", treat = "trt",
      base_cov = c("X_0", "W_0"),
      rx = "time_on_trt_relative",
      psi_test = "phreg",
      # psi_test = "lifereg",
      alpha = alpha,
      censor_time = "max_FU",
      autoswitch = TRUE,
      recensor = recensor,
      gridsearch = FALSE,
      root_finding = "bisection",
      boot = B > 0,
      n_boot = B
    )

    SE <- sd(log(RPS$hr_boots))
    list(
      HR = RPS$hr,
      SElogHR = SE,
      low = RPS$hr_CI[1],
      up = RPS$hr_CI[2],
      p = RPS$cox_pvalue
    )
  }
}
