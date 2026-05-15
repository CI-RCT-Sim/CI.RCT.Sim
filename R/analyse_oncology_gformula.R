#' Create Analyse Functions for g formula
#'
#' @param B Number of bootstrap samples
#' @param reps Number of repetitions for the g-formula simulation (per bootstrap sample)
#'
#' @return an analyse function that can be used in runSimulation
#' @export
#'
#' @importFrom survival tmerge coxph Surv
#' @importFrom stats confint glm predict sd model.matrix
#'
#' @examples
#' \donttest{
#' setting <- oncology_scenario()[1, ]
#'
#' dat <- generate_oncology(setting)
#'
#' analyse_oncology_gformula()(setting, dat)
#' }
analyse_oncology_gformula <- function(B = 20, reps = 1) {
  function(condition, dat, fixed_objects = NULL) {
    intervals_per_year <- 12
    data <- dat[order(dat$id), ]
    n <- dim(data)[1]

    IL <- 1 / intervals_per_year
    data_gform <- lapply(1:dim(data)[1], \(i){
      temp <- data[i, ]

      int_end <- floor(temp$event_time * intervals_per_year)
      D <- data.frame(id = i, time = 0:int_end, trt = temp$trt, ev = 0, prog = 0, switch = 0, X = temp$X_0, W = NA, L = NA, prog_now = 0, switch_now = 0, ever_prog = FALSE, ever_switch = FALSE) # time=interval
      D$start <- D$time * IL
      D$stop <- (D$time + 1) * IL

      # D$ev[temp$event_time > D$start & temp$event_time <= D$stop] <- 1
      if(temp$ev==1)  D$ev[int_end+1] <- 1

      D$prog[temp$prog_time <= D$stop] <- 1
      if (any(D$prog == 1)) {
        D$ever_prog <- TRUE
        D$prog_now[sum(D$prog == 0) + 1] <- 1
      }

      if (temp$switch == 1) {
        D$ever_switch <- TRUE
        D$switch <- D$prog
        D$switch_now <- D$prog_now
      }
      nD <- dim(D)[1]
      WW <- rep(as.numeric(temp[, paste("W", 0:(condition$k - 1), sep = "_")]), each = intervals_per_year)
      D$W <- WW[1:nD]
      LL <- rep(as.numeric(temp[, paste("L", 0:(condition$k - 1), sep = "_")]), each = intervals_per_year)
      D$L <- LL[1:nD]

      D$trt_actual <- D$trt + D$switch
      D$W_lag <- c(NA, D$W[-length(D$W)])
      D
    }) |>
      do.call(rbind, args=_)

    D_all <- D <- data_gform
    ids <- data$id
    hr <- rep(NA, B + 1)
    for (boot in 1:(B + 1)) {
      if (boot == 1) {
        D <- D_all
        dat_bs <- data
      } else {
        boot_ind <- sample(ids, n, replace = TRUE)
        D <- merge(D_all, data.frame(id=boot_ind), by="id", all.y = TRUE) |>
          sort_by(~id+time)
        dat_bs <- data[boot_ind, ] # data is orderd by id
      }

      mod_prog <- glm(prog ~ X + W + trt_actual + time, family = binomial, data = D[D$prog == 0 | D$prog_now == 1, ])
      mod_sw <- glm(switch ~ X + W, family = binomial, data = D[D$trt == 0 & D$prog_now == 1, ])
      mod_death <- glm(ev ~ X + W + trt + time + switch, family = binomial, data = D)
      mod_W <- lm(W ~ X + W_lag + time, data = D)

      DS <- NULL

      maxint <- round(condition$max_duration * intervals_per_year)

      switch <- rep(0, n)

      b_W <- coef(mod_W)
      sd_W <- summary(mod_W)$sigma
      b_death <- coef(mod_death)
      b_prog <- coef(mod_prog)
      b_switch <- coef(mod_sw)


      samp_binom <- function(X, b) {
        rbinom(dim(X)[1], size = 1, prob = plogis(X %*% b))
      }
      sim_with_switching <- FALSE

      ev_time <- rep(NA, n)

      active <- rep(TRUE, n)
      switch <- rep(0, n)
      prog <- rep(0, n)
      prog_now <- rep(0, n)
      trt <- dat_bs$trt
      i <- 1
      X <- X0 <- dat_bs$X_0
      W <- W0 <- dat_bs$W_0
      maxtime <- dat_bs$calendar_end_of_study - dat_bs$calendar_start_time
      maxint_i <- round(maxtime * intervals_per_year)
      event_final <- rep(0, n)
      TIME_EV <- lapply(1:reps, \(u){
        for (i in 1:maxint) {
          zeit <- i - 0.5 # i or i-0.5 to make it mid-interval

          M <- model.matrix(~ X + W + trt)
          trt_actual <- as.numeric(trt == 1 | switch == 1)
          M <- cbind(M, time = i - 1, switch = switch, trt_actual = trt_actual)
          if (i > 1) M[, "W"] <- M[, c(1:3, 5)] %*% b_W + rnorm(n, 0, sd_W)

          ev <- samp_binom(M[active, 1:6, drop = FALSE], b_death)
          set_ev <- ev == 1
          ev_time[active][set_ev] <- zeit
          event_final[active][set_ev] <- 1
          active[active][set_ev] <- FALSE
          if (!any(active)) break
          prog_now <- rep(0, n)
          prog_now[active] <- samp_binom(M[active, c(1:3, 7, 5), drop = FALSE], b_prog)
          prog[active] <- prog_now[active]
          if (sim_with_switching) {
            set_switch_possible <- active & trt == 0 & prog_now == 1
            if (any(set_switch_possible)) {
              switch[set_switch_possible] <- samp_binom(M[set_switch_possible, 1:3, drop = FALSE], b_switch)
            }
          }
          admin_cens <- active & maxint_i == i
          if (any(admin_cens)) {
            active[admin_cens] <- FALSE
            ev_time[admin_cens] <- zeit
          }
        }
        data.frame(TIME=ev_time, EV=event_final)
      }) |>
        do.call(rbind, args=_)

      X0 <- rep(X0, reps)
      W0 <- rep(W0, reps)
      trt <- rep(trt, reps)
      cox <- coxph(Surv(time = TIME_EV$TIME, event = TIME_EV$EV) ~ trt + X0 + W0)
      hr[boot] <- coef(cox)[1]
    }
    SE <- sd(hr[-1])
    p <- 2 * (1 - pnorm(abs(hr[1] / SE)))
    KI <- exp(hr[1] + c(-1, 1) * SE * qnorm(0.975))
    list(
      HR = exp(hr[1]),
      SElogHR = SE,
      low = KI[1],
      up = KI[2],
      p = p,
      N_pat = n,
      N_evt = sum(dat$ev)
    )
  }
}
