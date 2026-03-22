#' Generate Dataset that simulates an oncology trial
#'
#' @param condition condition row of Design dataset
#' @param fixed_objects fixed objects of Design dataset
#'
#' @details
#' Condition has to contain the following columns:
#'
#'   * k              # number of post baseline visits
#'   * recr_interval  # recruitment interval in years (e.g. 2 means that recruitment takes 2 years)
#'   * max_duration   # maximal duration of the trial in years (e.g. 7 means that the last patient recruited is followed for 5 years)
#'   * alpha          # type I error rate for the logrank test at the end of the trial (i.e. at max_duration)
#'   * power          # power for the logrank test at the end of the trial (i.e. at max_duration)
#'   * p_trt          # proportion of patients randomized to the treatment arm
#'   * w              # threshold for covariate W to determine the effect of W on the hazard
#'   * mu_W           # list of length 1 with two data.frames of dimension k x 2, one for control and one for treatment, with the mean values of W at each visit for control and treatment group
#'   * mu_L           # list of length 1 with two data.frames of dimension k x 2, one for control and one for treatment, with the mean values of L at each visit for control and treatment group
#'   * Sigma_W_L      # list of length 1 with one data.frame of dimension 2k x 2k with the covariance matrix of W and L at all visits (assuming the same covariance matrix for control and treatment group)
#'   * beta_prog      # list of length 1 with one data.frame of dimension 7 x k with the coefficients for the progression hazard at each visit for the variables Int, X, W, W>0, L, trt, switched
#'   * beta_switch    # list of length 1 with one data.frame of dimension 3 x k with the coefficients for the switching probability at each visit for the variables Int, X, W
#'   * beta_death     # list of length 1 with one data.frame of dimension 7 x k with the coefficients for the death hazard at each visit for the variables Int, X, W, W>0, L, trt, switched
#'   * beta_cens      # list of length 1 with one data.frame of dimension 7 x k with the coefficients for the censoring hazard at each visit for the variables Int, X, W, W>0, L, trt, switched
#'
#' @return
#' For generate_oncology: A data set with n rows and the columns .....
#'
#'
#' @export
#' @describeIn generate_oncology simulates a data set with n rows.
#'
#' @examples
#' Design <- oncology_scenario()
#' generate_oncology(Design[1, ])
generate_oncology <- function(condition, fixed_objects = list(allow_switch = TRUE, logHR_assumed = NULL, ev_soll = NULL)) {
  if (!(condition$k >= condition$max_duration + 1)) stop("k must be greater or equal max_duration+1")

  if (is.null(fixed_objects$logHR_assumed)) logHR_assumed <- as.numeric(condition$beta_death[[1]]["logHR_assumed"])
  #
  if (is.null(fixed_objects$ev_soll)) ev_soll <- ceiling(((qnorm(1 - condition$alpha / 2) + qnorm(condition$power)) / logHR_assumed)^2 / condition$p_trt / (1 - condition$p_trt))

  n_aim <- ceiling(ev_soll) * 2
  n <- rpois(1, n_aim)
  n0 <- rbinom(1, size = n, prob = condition$p_trt)
  n1 <- n - n0
  trt <- rep(c(0, 1), times = c(n0, n1))
  X <- matrix(rep(rnorm(n, 0, 1), k), ncol = k)
  W <- rbind(
    rmvnorm(n0, mean = condition$mu_W[[1]]$ctr, sigma = condition$Sigma_W_L[[1]]),
    rmvnorm(n1, mean = condition$mu_W[[1]]$trt, sigma = condition$Sigma_W_L[[1]])
  )
  Wgrw <- W > condition$w
  L <- rbind(
    rmvnorm(n0, mean = condition$mu_L[[1]]$ctr, sigma = condition$Sigma_W_L[[1]]),
    rmvnorm(n1, mean = condition$mu_L[[1]]$trt, sigma = condition$Sigma_W_L[[1]])
  )
  switchtime <- rep(Inf, n)

  prog_time <- rtime2(condition$beta_prog[[1]], k, X, W, Wgrw, L, trt, switchtime)

  # covariate values at progression (secondary baseline)
  index_sec_BL <- floor(prog_time) + 1
  index_sec_BL[index_sec_BL > k] <- k

  X_2BL <- W_2BL <- L_2BL <- TRT_2BL <- rep(NA, n)
  for (j in 1:n) {
    X_2BL[j] <- X[j, index_sec_BL[j]]
    W_2BL[j] <- W[j, index_sec_BL[j]]
    L_2BL[j] <- L[j, index_sec_BL[j]]
  }

  h_0 <- condition$beta_switch[[1]]["Int"]
  h_X <- condition$beta_switch[[1]]["X"]
  h_W <- condition$beta_switch[[1]]["W"]

  if (fixed_objects$allow_switch) {
    switch_prob <- ifelse(trt == 0, expit(h_0 + X_2BL * h_X + W_2BL * h_W), 0)
  } else {
    switch_prob <- rep(0, n)
  }
  switch <- rbinom(n, size = 1, prob = switch_prob) == 1
  switchtime[switch] <- prog_time[switch]

  TRT_2BL <- trt
  TRT_2BL[switch] <- 1

  ind <- 0:(k - 1) # used for variable names to indicate visits. Visit 0 is baseline, then visits are performed every year. So ind is the time of the visit in years.

  # event_time_uncensored<-rtime(n,TRT,SWITCHED,X,W,L,b0,b_trt,b_sw,b_x,b_W,b_L)
  event_time_uncensored <- rtime2(condition$beta_death[[1]][1:7], k, X, W, Wgrw, L, trt, switchtime)
  if (condition$beta_cens[[1]]["Int"] > -Inf) {
    random_cens_time <- rtime2(condition$beta_cens[[1]], k, X, W, Wgrw, L, trt, switchtime)
  } else {
    # this is to save time, the rtime2 function would return Inf, too
    random_cens_time <- rep(Inf, n)
  }

  rcens <- random_cens_time < event_time_uncensored
  # mean(random_cens_time>1)
  # mean(rcens)
  event_time_rc <- ifelse(rcens, random_cens_time, event_time_uncensored)
  event <- as.numeric(!rcens)
  start <- runif(n, 0, condition$recr_interval)
  cal <- start + event_time_rc # uncensored #calendar times


  colnames(X) <- paste("X", ind, sep = "_")
  colnames(W) <- paste("W", ind, sep = "_")
  colnames(L) <- paste("L", ind, sep = "_")


  temp <- data.frame(
    id = 1:n,
    trt,
    # TRT,
    X, W, L,
    # PD,switch,SWITCHED,index_sec_BL,TRT_2BL,
    X_2BL, W_2BL, L_2BL,
    cal,
    event_time_rc,
    start,
    event_time = NA,
    ev = event,
    prog_time,
    prog_ev = NA,
    switch
  )

  # administrative censoring
  temp <- temp[order(temp$cal), ]
  # temp$ev<-1
  temp$cum_ev <- cumsum(temp$ev)
  index_end <- sum(temp$cum_ev <= ev_soll)
  cal_end <- temp$cal[index_end]

  temp$admin_cens <- temp$cal > cal_end
  temp$ev[temp$admin_cens] <- 0 # as.numeric(temp$cal<=cal_end)
  temp$event_time <- ifelse(temp$admin_cens, cal_end - temp$start, temp$event_time_rc)
  temp$prog_ev <- as.numeric(temp$prog_time < temp$event_time)
  temp$prog_time[!temp$prog_ev] <- temp$event_time[!temp$prog_ev]
  temp$calendar_start_time <- temp$start
  temp$calendar_end_of_study <- cal_end

  temp$cal <- NULL
  temp$event_time_rc <- NULL

  temp$start <- NULL
  temp$cum_ev <- NULL
  temp$cal2 <- NULL

  temp$switch <- as.numeric(temp$switch & temp$prog_ev == 1)

  # remove covariate values at unobserved timepoints
  ind <- 0:(k - 1)
  i <- 1
  for (i in 1:dim(temp)[1]) {
    unobs <- ind > temp$event_time[i]
    rm_X <- paste(rep(c("X", "W", "L"), each = sum(unobs)), ind[unobs], sep = "_")
    temp[i, rm_X] <- NA
  }
  temp$random_cens <- temp$ev == 0 & !temp$admin_cens
  temp$ev_soll <- ev_soll
  temp$ev_obs <- sum(temp$ev)

  temp
}

#' Create an empty assumptions data.frame for generate_oncology
#'
#' @param print print code to generate parameter set?
#'
#' @return For oncology_scenario: a design tibble with default values invisibly
#'
#' @details oncology_scenario generates a default design `data.frame`
#'   for use with generate_oncology. If print is `TRUE` code to produce
#'   the template is also printed for copying, pasting and editing by the user.
#'   (This is the default when run in an interactive session.)
#'
#' @export
#' @describeIn generate_oncology generate default design tibble
#'
#' @examples
#' Design <- oncology_scenario()[1, ]
#' Design
oncology_scenario <- function(print = interactive()) {
  skel <- "params_scenarios_grid(
  k             = 10,     # Number of visits post baseline
  recr_interval = 2,      # Recruitment interval in years (e.g. 2 means that recruitment takes 2 years)
  max_duration  = 7,      # Maximal duration of the trial in years (e.g. 7 means that the last patient recruited is followed for 5 years)
  alpha         = 0.05,   # Type I error rate for the logrank test at the end of the trial (i.e. at max_duration)
  power         = 0.8,    # Power for the logrank test at the end of the trial (i.e. at max_duration)
  p_trt         = 0.5,    # Proportion of patients randomized to the treatment arm
  w             = 0,      # Threshold for covariate W to determine the effect of W on the hazard
  mu_W          =
  mu_L          =
  Sigma_W_L
  beta_prog
  beta_switch
  beta_death
  beta_cens


  ) |>
  merge(data.frame(hyp=c(1,0)), by=NULL)
"


  if (print) {
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
#' @describeIn generate_oncology  calculate true summary statistics for ...
#'
#' @examples
#' oncology_scenario_set_truevalues(oncology_scenario())
oncology_scenario_set_truevalues <- function(Design, cutoff_stats = 10, fixed_objects = NULL) {
  # true_summary_statistics_diabetes_rescue_rowwise <- function(condition, cutoff_stats) {
  #   res <- data.frame(
  #     eff_true <- condition$delta / 2 * (1 - exp(-condition$lambda * condition$k))
  #   )
  #   res
  # }
  #
  # Design <- Design |>
  #   split(1:nrow(Design)) |>
  #   mapply(FUN = true_summary_statistics_diabetes_rescue_rowwise, cutoff_stats = cutoff_stats, SIMPLIFY = FALSE)
  #
  # Design <- do.call(rbind, Design)
  Design$eff_true <- Design$delta / 2 * (1 - exp(-Design$lambda * Design$k))

  # specifying parameters for sample size calculation
  alpha <- 0.05
  power <- 0.8

  if (!"nfix" %in% colnames(Design)) {
    Design$n <- 2 * round(((qnorm(1 - alpha / 2) + qnorm(power))^2) *
      Design$sd_bl^2 * (1 - Design$rho^2) * 2 /
      (Design$eff_true^2))
  } else {
    Design$n <- Design$nfix
  }

  Design$eff_true <- ifelse(Design$hyp == 1, Design$eff_true, 0)
  Design
}
