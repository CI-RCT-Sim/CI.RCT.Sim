#' Create Analyse Functions for
#'
#' @param use_censoring_IPW logical, Default is FALSE
#' @param trunc_weights numeric
#' @param requ_n_cens numeric
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
#' analyse_oncology_ipw2()(setting, dat)
analyse_oncology_ipw2 <- function(use_censoring_IPW = FALSE,
                                  trunc_weights = 5,
                                  requ_n_cens = 5) {
  function(condition, dat, fixed_objects = NULL) {
    prep_data_RPSFTM_fun <- function(data) {
      data$time_on_trt <- 0
      data$time_on_trt[data$trt == 1] <- data$event_time[data$trt == 1]
      set_sw <- data$trt == 0 & data$switch == 1
      data$time_on_trt[set_sw] <- data$event_time[set_sw] - data$prog_time[set_sw]
      data$time_on_trt_relative <- data$time_on_trt / data$event_time
      data$max_FU <- data$calendar_end_of_study - data$calendar_start_time
      set_admin_cens <- data$max_FU < data$event_time
      data$max_FU[set_admin_cens] <- data$event_time[set_admin_cens]
      data
    }

    set <- dat$trt == 0 & dat$prog_ev == 1
    wmod <- glm(
      switch ~ X_2BL + W_2BL,
      family = binomial,
      data = dat,
      subset = set
    )
    pred <- predict(wmod, type = "response")
    pred_a <- ifelse(dat$switch[set] == 1, pred, 1 - pred)
    dat$w <- 1
    dat$w[set] <- 1 / pred_a


    #only required columns
    data_k <- dat[, c("id",
                      "trt",
                      "X_0",
                      "W_0",
                      "event_time",
                      "ev",
                      "prog_time",
                      "switch",
                      "w")]
    sdat <- tmerge(
      data1 = data_k,
      data2 = data_k,
      id = id,
      tstop = event_time
    )
    sdat <- tmerge(
      data1 = sdat,
      data2 = data_k,
      id = id,
      death_event = event(event_time, ev)
    )
    sdat <- tmerge(
      data1 = sdat,
      data2 = data_k,
      id = id,
      PD = tdc(prog_time)
    )

    sdat$w[sdat$PD == 0 | sdat$trt == 1] <- 1
    sdat$w[sdat$w > trunc_weights] <- trunc_weights

    if (use_censoring_IPW & sum(dat$random_cens) >= requ_n_cens) {
      long <- reshape(
        data[, c(
          "id",
          "trt",
          "X_0",
          paste("W", 0:9, sep = "_"),
          "event_time",
          "ev",
          "prog_time",
          "switch",
          "admin_cens",
          "random_cens"
        )],
        direction = "long",
        varying = paste("W", 0:9, sep = "_"),
        sep = "_"
      )
      long <- long[!is.na(long$W), ]
      long <- long[order(long$id, long$time), ] #important to sort by id, so it is in this order later

      long$start <- long$time
      long$stop <- long$time + 1
      long$random_cens_event <- 0
      ids <- unique(data$id)

      for (id in ids) {
        set <- long$id == id
        temp <- long[set, ]
        m <- dim(temp)[1]
        temp$stop[m] <- temp$event_time[m]
        if (m > 1)
          temp$ev[1:(m - 1)] <- 0
        temp$random_cens_event[m] <- as.numeric(temp$random_cens[m])
        long[set, ] <- temp
      }


      tdmod <- coxph(
        Surv(
          time = start,
          time2 = stop,
          event = random_cens_event
        ) ~ trt * (X_0 + W),
        id = id,
        data = long
      ) #there is a warning, if only one group contains random censoring, but this does not matter

      fit <- survfit(tdmod, id = id, newdata = long)

      remainprob <- data.frame(
        id = rep(as.numeric(names(fit$strata)), times = fit$strata),
        time = fit$time,
        surv = fit$surv
      )
      remainprob$IPW <- 1 / remainprob$surv

      rund <- 12
      sdat$tstart <- round(sdat$tstart, rund)
      sdat$tstop <- round(sdat$tstop, rund)
      remainprob$time <- round(remainprob$time, rund)

      sdat <- tmerge(
        data1 = sdat,
        data2 = remainprob,
        id = id,
        censIPW = tdc(time, IPW)
      )

      sdat$censIPW[is.na(sdat$censIPW)] <- 1 #for time 0 to first cens/event time there is now entry in the fitted curve, need to fill it here
      sdat$censIPW[sdat$censIPW > trunc_weights] <- trunc_weights #to avoid extreme weights. These may occur if there are few censoring events
      sdat$w <- sdat$w * sdat$censIPW
    }
    sdat$w[sdat$w > trunc_weights] <- trunc_weights

    remo <- sdat$trt == 0 & sdat$PD == 1 & sdat$switch == 1
    dat_a <- sdat[!remo, ]

    mod <- coxph(
      Surv(
        time = tstart,
        time2 = tstop,
        event = death_event
      ) ~ trt + X_0 + W_0,
      data = dat_a,
      weights = w,
      robust = TRUE,
      id = id
    )

    HR <- exp(coef(mod)[["trt"]])
    CI <- exp(confint(mod)[1, ])
    smr <- summary(mod)
    p <- smr$coef[1, "Pr(>|z|)"]
    SE <- smr$coef[1, "robust se"]
    list(
      HR = HR,
      SElogHR = SE,
      low = CI[["2.5 %"]],
      up = CI[["97.5 %"]],
      p = p
    )
  }
}
