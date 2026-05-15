#' Create Analyse Functions for
#'
#' @param method character, method to adjust for treatment switching, either "RPSFTM" or "TSE". Default is "RPSFTM".
#' @param recensor logical, whether to recensor the data after adjusting for treatment switching. Default is TRUE.
#' @param alpha numeric, significance level for confidence intervals. Default is 0.05.
#' @param B integer, number of bootstrap samples for estimating standard errors. Default is 100.
#' @param trunc_weights numeric, maximum value for inverse probability weights to avoid extreme weights. Default is 5.
#' @param show_progress logical, whether to print progress during bootstrapping. Default is FALSE.
#' @param use_censoring_IPW logical, whether to use inverse probability weighting to adjust for censoring in the final Cox model. Default is FALSE.
#' @param requ_n_cens numeric, minimum number of censoring events required to use inverse probability weighting for censoring. Default is 5.
#'
#' @return an analyse function that can be used in runSimulation
#' @export
#'
#' @importFrom trtswitch rpsftm tsesimp
#'
#' @examples
#' setting <- oncology_scenario()[1, ]
#' dat <- generate_oncology(setting)
#' analyse_oncology_mixed(method = "RPSFTM")(setting, dat)
#' analyse_oncology_mixed(method = "TSE")(setting, dat)
analyse_oncology_mixed <- function(method = c("RPSFTM", "TSE")[1],
                                   recensor = TRUE,
                                   alpha = 0.05,
                                   B = 100,
                                   trunc_weights = 5,
                                   show_progress = FALSE,
                                   use_censoring_IPW = FALSE,
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
    get_weights<-function(data_k,trunc_weights=Inf) {

      sdat<- tmerge(data1=data_k, data2=data_k, id=id, tstop=event_time)
      sdat<- tmerge(data1=sdat, data2=data_k, id=id, death_event = event(event_time,ev))

      long<-reshape(data[,c("id","trt","X_0",paste("W",0:9,sep="_"),"event_time","ev","prog_time","switch","admin_cens","random_cens")],
                    direction = "long", varying = paste("W",0:9,sep="_"), sep="_")
      long<-long[!is.na(long$W),]
      long<-long[order(long$id,long$time),] #important to sort by id, so it is in this order later

      long$start<-long$time
      long$stop<-long$time+1
      long$random_cens_event<-0
      ids<-unique(data$id)

      for(id in ids) {
        set<-long$id==id
        temp<-long[set,]
        m<-dim(temp)[1]
        temp$stop[m]<-temp$event_time[m]
        if(m>1) temp$ev[1:(m-1)]<-0
        temp$random_cens_event[m]<-as.numeric(temp$random_cens[m])
        long[set,]<-temp
      }

      table(long$random_cens_event,long$trt)

      tdmod<-coxph(Surv(time=start,time2=stop,event=random_cens_event)~trt*(X_0+W),id=id,data=long)
      fit<-survfit(tdmod,id=id,newdata=long)

      remainprob<-data.frame(id=rep(as.numeric(names(fit$strata)),times=fit$strata),time=fit$time,surv=fit$surv)
      remainprob$IPW<-1/remainprob$surv

      rund<-12
      sdat$tstart<-round(sdat$tstart,rund)
      sdat$tstop<-round(sdat$tstop,rund)
      remainprob$time<-round(remainprob$time,rund)

      sdat<-tmerge(data1=sdat, data2=remainprob, id=id, censIPW = tdc(time, IPW))

      sdat$censIPW[is.na(sdat$censIPW)]<-1 #for time 0 to first cens/event time there is now entry in the fitted curve, need to fill it here
      sdat$censIPW[sdat$censIPW>trunc_weights]<-trunc_weights #to avoid extreme weights. These may occur if there are few censoring events
      sdat$w<-sdat$censIPW

      sdat$w[sdat$w>trunc_weights]<-trunc_weights

      sdat
    }
    data <- dat[order(dat$id), ]
    n <- dim(data)[1]

    if (B < 1)
      B <- 1 #the first sample is always the original data

    ids <- data$id

    hr <- rep(NA, B + 1) #actuall, it will be filled with log(HR)

    for (boot in 1:(B + 1)) {
      if (show_progress) {
        print(boot)
        utils::flush.console()
      }
      if (boot == 1) {
        dat_bs <- data
      } else {
        boot_ind <- sample(ids, n, replace = TRUE)
        dat_bs <- data[boot_ind, ] #data is orderd by id
        dat_bs$id <- 1:dim(dat_bs)[1]
      }

      data_rpsftm <- prep_data_RPSFTM_fun(dat_bs)

      ###########
      if (method == "RPSFTM") {
        analysis <- rpsftm(
          data = data_rpsftm,
          id = "id",
          time = "event_time",
          event = "ev",
          treat = "trt",
          base_cov = c("X_0", "W_0"),
          ##
          rx = "time_on_trt_relative",
          psi_test = "phreg",
          alpha = alpha,
          censor_time = "max_FU",
          autoswitch = TRUE,
          recensor = recensor,
          gridsearch = FALSE,
          root_finding = "bisection",
          boot = FALSE
        )
      }
      if (method == "TSE") {
        analysis <- tsesimp(
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
          boot = FALSE
        )
      }

      data_k <- analysis$data_outcome
      names(data_k)[names(data_k) == "t_star"] <- "event_time"
      names(data_k)[names(data_k) == "d_star"] <- "ev"

      if (use_censoring_IPW & sum(data$random_cens) >= requ_n_cens) {
        sdat <- get_weights(data_k, trunc_weights = trunc_weights)
        cox <- coxph(
          Surv(
            time = tstart,
            time2 = tstop,
            event = death_event
          ) ~ trt + X_0 + W_0,
          data = sdat,
          weights = w,
          robust = TRUE,
          id = id
        )
      } else {
        cox <- coxph(Surv(time = event_time, event = ev) ~ trt + X_0 + W_0, data =
                       data_k)
      }
      hr[boot] <- coef(cox)[1]
    }

    if (B > 0) {
      SE <- sd(hr[-1])
      p <- 2 * (1 - pnorm(abs(hr[1] / SE)))
      KI <- exp(hr[1] + c(-1, 1) * SE * qnorm(0.975))
      out <- list(
        HR = exp(hr[1]),
        SElogHR = SE,
        low = KI[1],
        up = KI[2],
        p = p
      )
    } else {
      out <- list(
        HR = exp(hr[1]),
        SElogHR = NA,
        low = NA,
        up = NA,
        p = NA
      )
    }
    out
  }
}
