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
#'   * aimed_for_n_per_required_event # This governs the recruitment rate.
#'   * alpha          # type I error rate for the logrank test at the end of the trial (i.e. at max_duration)
#'   * power          # power for the logrank test at the end of the trial (i.e. at max_duration)
#'   * p_trt          # proportion of patients randomized to the treatment arm
#'   * w              # threshold for covariate W to determine the effect of W on the hazard
#'   * mu_W           # list with list entries containing mean values of W at each visit for control and treatment group in k dimensional vectors
#'   * mu_L           # list with list entries containing mean values of L at each visit for control and treatment group in k dimensional vectors
#'   * Sigma_W_L      # list with matrix of dimension k x k with the covariance matrix of W and L at all visits (assuming the same covariance matrix for control and treatment group)
#'   * beta_prog      # list with vectors of dimension 7 with the coefficients for the progression hazard at each visit for the variables Int, X, W, W>0, L, trt, switched
#'   * beta_switch    # list with vectors of dimension 7 with the coefficients for the switching probability at each visit for the variables Int, X, W, W>0, L, trt, switched
#'   * beta_cens      # list with vectors of dimension 7 with the coefficients for the censoring hazard at each visit for the variables Int, X, W, W>0, L, trt, switched
#'   * beta_death     # list with vectors of dimension 8 with the coefficients for the death hazard at each visit for the variables Int, X, W, W>0, L, trt, switched, logHR_assumed
#'
#' @return
#' For generate_oncology: A data set with n rows and the columns .....
#'
#' @importFrom stats plogis rbinom rpois runif rnorm
#' @importFrom SimDesign rmvnorm
#'
#' @export
#' @describeIn generate_oncology simulates a data set with n rows.
#'
#' @examples
#' Design <- oncology_scenario()
#' generate_oncology(Design[1, ])
generate_oncology <- function(condition, fixed_objects = list(allow_switch = TRUE, logHR_assumed = NULL, ev_soll = NULL, allow_random_cens = TRUE)) {
  if (!(condition$k >= (condition$max_duration + 1))) stop("k must be greater or equal max_duration+1")

  names(condition$beta_prog[[1]]) <-
    names(condition$beta_switch[[1]]) <-
    c("Int", "X", "W", "Wgrw", "L", "trt", "switched")


  names(condition$beta_cens[[1]]) <-
    c("Int", "X", "W", "Wgrw", "L", "trt", "switched", "control_only_random_cens")

  names(condition$beta_death[[1]]) <-
    c("Int", "X", "W", "Wgrw", "L", "trt", "switched", "logHR_assumed")

  rtime2 <- function(b, k, X, W, Wgrw, L, trt, switchtime) {
    n <- dim(X)[1]
    time <- rep(NA, n)
    u <- -log(runif(n))
    for (i in 1:n) {
      temp <- data.frame(time = 0:(k - 1), Int = 1, X = X[i, ], W = W[i, ], Wgrw = Wgrw[i, ], L = L[i, ], trt = trt[i])
      if (switchtime[i] < Inf) {
        ind <- sum(temp$time <= switchtime[i])
        newline <- temp[ind, ]
        newline$time <- switchtime[i]
        temp <- rbind(temp, newline)
      }
      temp <- temp[order(temp$time), ]
      temp$switched <- as.numeric(temp$time >= switchtime[i])


      temp$haz <- exp(as.matrix(temp[, -1]) %*% b)
      ntime <- dim(temp)[1]
      temp$cumhaz <- c(0, cumsum(diff(temp$time) * temp$haz[-ntime]))

      a <- sum(temp$cumhaz <= u[i])
      if (a >= ntime) {
        time[i] <- temp$time[ntime] + (u[i] - temp$cumhaz[ntime]) / temp$haz[ntime]
      } else {
        temp$cumhaz[a]
        s <- (u[i] - temp$cumhaz[a]) / (temp$cumhaz[a + 1] - temp$cumhaz[a])
        time[i] <- temp$time[a] + s * (temp$time[a + 1] - temp$time[a])
      }
    }
    time
  }

  ifelse(is.null(fixed_objects$logHR_assumed),
    logHR_assumed <- as.numeric(condition$beta_death[[1]]["logHR_assumed"]),
    logHR_assumed <- fixed_objects$logHR_assumed
  )

  ifelse(is.null(fixed_objects$ev_soll),
    ev_soll <- ceiling(((qnorm(1 - condition$alpha / 2) + qnorm(condition$power)) / logHR_assumed)^2 / condition$p_trt / (1 - condition$p_trt)),
    ev_soll <- fixed_objects$ev_soll
  )

  n_aim <- ceiling(ev_soll) * condition$aimed_for_n_per_required_event
  n <- rpois(1, n_aim)
  n0 <- rbinom(1, size = n, prob = condition$p_trt)
  n1 <- n - n0
  trt <- rep(c(0, 1), times = c(n0, n1))
  X <- matrix(rep(rnorm(n, 0, 1), condition$k), ncol = condition$k)

  W <- rbind(
    rmvnorm(n0, mean = condition$mu_W[[1]]$ctr, sigma = as.matrix(condition$Sigma_W_L[[1]])),
    rmvnorm(n1, mean = condition$mu_W[[1]]$trt, sigma = condition$Sigma_W_L[[1]])
  )
  Wgrw <- W > condition$w
  L <- rbind(
    rmvnorm(n0, mean = condition$mu_L[[1]]$ctr, sigma = condition$Sigma_W_L[[1]]),
    rmvnorm(n1, mean = condition$mu_L[[1]]$trt, sigma = condition$Sigma_W_L[[1]])
  )
  switchtime <- rep(Inf, n)

  prog_time <- rtime2(condition$beta_prog[[1]], condition$k, X, W, Wgrw, L, trt, switchtime)

  # covariate values at progression (secondary baseline)
  index_sec_BL <- floor(prog_time) + 1
  index_sec_BL[index_sec_BL > condition$k] <- condition$k

  X_2BL <- W_2BL <- L_2BL <- TRT_2BL <- Wgrw_2BL <- rep(NA, n)
  for (j in 1:n) {
    X_2BL[j] <- X[j, index_sec_BL[j]]
    W_2BL[j] <- W[j, index_sec_BL[j]]
    L_2BL[j] <- L[j, index_sec_BL[j]]
    Wgrw_2BL[j] <- Wgrw[j, index_sec_BL[j]]
  }

  h_0 <- condition$beta_switch[[1]]["Int"]
  h_X <- condition$beta_switch[[1]]["X"]
  h_W <- condition$beta_switch[[1]]["W"]
  h_Wgrw <- condition$beta_switch[[1]]["Wgrw"]
  h_L <- condition$beta_switch[[1]]["L"]

  if (fixed_objects$allow_switch) {
    switch_prob <-
      ifelse(trt == 0,
        plogis(h_0 + X_2BL * h_X + W_2BL * h_W + Wgrw_2BL * h_Wgrw + L_2BL * h_L),
        0
      )
  } else {
    switch_prob <- rep(0, n)
  }
  switch <- rbinom(n, size = 1, prob = switch_prob) == 1
  switchtime[switch] <- prog_time[switch]

  TRT_2BL <- trt
  TRT_2BL[switch] <- 1

  ind <- 0:(condition$k - 1) # used for variable names to indicate visits. Visit 0 is baseline, then visits are performed every year. So ind is the time of the visit in years.

  # event_time_uncensored<-rtime(n,TRT,SWITCHED,X,W,L,b0,b_trt,b_sw,b_x,b_W,b_L)
  event_time_uncensored <- rtime2(condition$beta_death[[1]][1:7], condition$k, X, W, Wgrw, L, trt, switchtime)
  if (condition$beta_cens[[1]]["Int"] > -Inf) {
    random_cens_time <- rtime2(condition$beta_cens[[1]][1:7], condition$k, X, W, Wgrw, L, trt, switchtime)
  } else {
    # this is to save time, the rtime2 function would return Inf, too
    random_cens_time <- rep(Inf, n)
  }

  # rcens <- random_cens_time < event_time_uncensored
  ## mean(random_cens_time>1)
  ## mean(rcens)
  # event_time_rc <- ifelse(rcens, random_cens_time, event_time_uncensored)
  # event <- as.numeric(!rcens)
  if (fixed_objects$allow_random_cens) {
    rcens <- random_cens_time < event_time_uncensored & (trt * condition$beta_cens[[1]]["control_only_random_cens"]) == 0
    # mean(random_cens_time>1)
    # mean(rcens)
    event_time_rc <- ifelse(rcens, random_cens_time, event_time_uncensored)
    event <- as.numeric(!rcens)
  } else {
    event_time_rc <- event_time_uncensored
    event <- rep(1, n)
  }


  start <- runif(n, 0, condition$recr_interval)
  cal <- start + event_time_rc # uncensored with respect to administrative censoring #calendar times


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
  # apply maximum study duration
  if (cal_end > condition$max_duration) cal_end <- condition$max_duration


  temp$admin_cens <- temp$cal > cal_end
  temp$ev[temp$admin_cens] <- 0 # as.numeric(temp$cal<=cal_end)
  temp$event_time <- ifelse(temp$admin_cens, cal_end - temp$start, temp$event_time_rc)

  # remove patients who were recruited after cal_end
  temp <- temp[temp$event_time > 0, ]

  #reset id, to have consecutive ID numbers
  temp$id<-1:dim(temp)[1]

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
  ind <- 0:(condition$k - 1)
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
#' @return For oncology_scenario: a design tibble with default values invisibly
#'
#' @details oncology_scenario generates a default design `data.frame`
#'   for use with generate_oncology.
#'
#' @export
#' @describeIn generate_oncology generate default design tibble
#'
#' @examples
#' Design <- oncology_scenario()[1, ]
#' Design
oncology_scenario <- function() {

  params_scenarios_grid <- function(...){
    params_args <- list(...)
    params_ref <- purrr::map(params_args, \(x){
      x[1]
    }) |>
      tibble::as_tibble()

    params_other <- purrr::imap(params_args, \(x, i){
      if(inherits(x, "list")){
        purrr::map(x[-1], \(y){
          tmp <- params_ref
          tmp[,i] <- list(list(y))
          tmp
        }) |>
          purrr::list_rbind()
      } else {
        purrr::map(x[-1], \(y){
          tmp <- params_ref
          tmp[,i] <- y
          tmp
        }) |>
          purrr::list_rbind()
      }
    }) |>
      purrr::list_rbind()

    rbind(params_ref, params_other)
  }




  AR_1_fun<-function(rho,k) {
    i<-rep(1:k,k)
    j<-rep(1:k,each=k)
    d<-abs(i-j)
    corr<-rho^d
    matrix(corr,k,k)
  }
  #AR_1_fun(0.9,10)

  exch_fun<-function(rho,k) {
    M<-matrix(rho,k,k)
    diag(M)<-1
    M
  }

  toeplitz_fun<-function(decr=0.1,k) {
    korr_vek<-seq(from=1, by=-decr, length.out=k)
    korr_vek[korr_vek<0]<-0
    i<-rep(1:k,k)
    j<-rep(1:k,each=k)
    d<-abs(i-j)+1
    corr<- korr_vek[d]
    matrix(corr,k,k)
  }
  #toeplitz_fun(0.1,10)


  k<-10

  label_fun<-function(x,labels) {
    for(i in 1:length(x)) names(x[[i]])<-labels
    x
  }

  mu_W<-list(
    list(trt=rep(0,k),ctr=rep(0,k)),
    list(trt=rep(0.5,k),ctr=c(0.5,0,-0.5,rep(-1,k-3))),
    list(trt=c(1,0.5,0,rep(-1,k-3)),ctr=c(1,0.5,0,rep(-1,k-3)))
  )
  mu_L<-mu_W
  Sigma_W_L<-list(exch_fun(0.5,k),toeplitz_fun(0.1,k))

  beta_lab<-c("Int","X","W","Wgrw","L","trt","switched")
  beta_lab2<-c("Int","X","W","Wgrw","L","trt","switched","logHR_assumed")
  beta_lab3<-c("Int","X","W","Wgrw","L","trt","switched","switching_in_control_only")


  beta_prog<-list(
    #    Int,             X,          W,    W>0, L,       trt ,      switched
    c( log( log(2)/0.5), log(0.5), log(0.5),  0, 0,        log(0.5),  0 ),
    c( log( log(2)/1)  , log(0.5), log(0.5),  0, 0,        log(0.5),  0 ),
    c( log( log(2)/0.5), log(0.5), 0       ,  0, 0,        log(0.5),  0 ),
    c( log( log(2)/0.5), log(0.5), log(0.5),  0, 0,        log(0.75), 0 ),
    c( log( log(2)/0.5), log(0.5), log(0.5),  0, 0,        0,         0 ),
    c( log( log(2)/0.5), log(0.5), log(0.5),  0, log(0.5), log(0.5),  0 )
  ) |> label_fun(beta_lab)

  beta_switch<-list(
    #   Int,          X,       W,      W>0,                      L   ,trt,  ,switched
    c(log(0.5/0.5),log(1.5),log(1.5),0                         ,0       ,0,0),
    c(log(0.9/0.1),log(1.5),log(1.5),0                         ,0       ,0,0),
    c(log(0.5/0.5),log(1.5),0       ,0                         ,0       ,0,0),
    c(log(0.5/0.5),log(1.5),log(1.5),log(0.9/0.1)*sqrt(pi/2)-log(1.5),0       ,0,0),
    c(log(0.5/0.5),log(1.5),log(1.5),0                         ,log(1.5),0,0)
  ) |> label_fun(beta_lab)

  HR_h<-0.5 #HR assumed high effect

  beta_death_high_effect<-list(
    #     Int,          X,          W,    W>0,  L,          trt,        switched,   HR_assumed
    c( log( log(2)/4), log(0.5), log(0.5),  0, 0       ,  log(0.5) ,    log(0.5) , log(HR_h) ),
    c( log( log(2)/2), log(0.5), log(0.5),  0, 0       ,  log(0.5) ,    log(0.5) , log(HR_h) ),
    c( log( log(2)/4), log(0.5), 0       ,  0, 0       ,  log(0.5) ,    log(0.5) , log(HR_h) ),
    c( log( log(2)/4), log(0.5), log(0.5),  0, log(0.5),  log(0.5) ,    log(0.5) , log(HR_h) ),
    c( log( log(2)/4), log(0.5), log(0.5),  0, 0       ,  log(0.5) ,    log(0.75), log(HR_h) ),
    c( log( log(2)/4), log(0.5), log(0.5),  0, 0       ,  log(0.5) ,    0        , log(HR_h) )
  ) |> label_fun(beta_lab2)


  HR_l<-0.75 #HR assumed low effect
  beta_death_low_effect<-list(
    #     Int,          X,          W,    W>0,  L,          trt,        switched   ,   HR_assumed
    c( log( log(2)/4), log(0.5), log(0.5),  0, 0       ,  log(0.75) ,    log(0.75) , log(HR_l)  ),
    c( log( log(2)/2), log(0.5), log(0.5),  0, 0       ,  log(0.75) ,    log(0.75) , log(HR_l)  ),
    c( log( log(2)/4), log(0.5), 0       ,  0, 0       ,  log(0.75) ,    log(0.75) , log(HR_l)  ),
    c( log( log(2)/4), log(0.5), log(0.5),  0, log(0.5),  log(0.75) ,    log(0.75) , log(HR_l)  ),
    c( log( log(2)/4), log(0.5), log(0.5),  0, 0       ,  log(0.75) ,    0         , log(HR_l)  )
  ) |> label_fun(beta_lab2)



  beta_death_H0_high_effect<-list(
    #    Int,              X,       W,    W>0,   L,       trt, switched   , HR_assumed
    c( log( log(2)/4), log(0.5), log(0.5),  0, 0       ,  0 ,    0 , log(HR_h)  ),
    c( log( log(2)/2), log(0.5), log(0.5),  0, 0       ,  0 ,    0 , log(HR_h)  ),
    c( log( log(2)/4), log(0.5), 0       ,  0, 0       ,  0 ,    0 , log(HR_h)  ),
    c( log( log(2)/4), log(0.5), log(0.5),  0, log(0.5),  0 ,    0 , log(HR_h)  )
  ) |> label_fun(beta_lab2)



  beta_death_H0_low_effect<-list(
    #    Int,              X,       W,     W>0,  L,       trt, switched   , HR_assumed
    c( log( log(2)/4), log(0.5), log(0.5),  0, 0       ,  0 ,    0 , log(HR_l)  ),
    c( log( log(2)/2), log(0.5), log(0.5),  0, 0       ,  0 ,    0 , log(HR_l)  ),
    c( log( log(2)/4), log(0.5), 0       ,  0, 0       ,  0 ,    0 , log(HR_l)  ),
    c( log( log(2)/4), log(0.5), log(0.5),  0, log(0.5),  0 ,    0 , log(HR_l)  )
  )  |> label_fun(beta_lab2)


  beta_cens<-list(
    #    Int,                X,       W,       W>0, L,        trt, switched, switching in control only
    c(log(-log(1-0.025))  ,0       , 0       ,  0,  0       ,    0,  0,		0),
    c(-Inf                ,0       , 0       ,  0,  0       ,    0,  0,		0),
    c(log(-log(1-0.025))  ,log(0.5), 0       ,  0,  0       ,    0,  0,		0),
    c(log(-log(1-0.025))  ,log(0.5), log(0.5),  0,  0       ,    0,  0,		0),
    c(log(-log(1-0.025))  ,log(0.5), log(0.5),  0,  log(0.5),    0,  0,		0),
    #extra
    c(log(-log(1-0.1))  ,log(0.5), log(0.5),  0,  0       ,    0,  0,			0),
    c(log(-log(1-0.1))  ,log(0.5), log(0.5),  0,  0       ,    0,  0,			1)


  ) |> label_fun(beta_lab3)



  constant<-data.frame(
    k=k,
    recr_interval=2,
    max_duration=7,
    alpha=0.05,
    power=0.8,
    p_trt=0.5,
    w=0,
    aimed_for_n_per_required_event=3,
    version=1
  )



  scen_trt_sw_H1_high<-cbind(params_scenarios_grid(mu_W=mu_W,mu_L=mu_L,Sigma_W_L=Sigma_W_L,beta_prog=beta_prog,beta_switch=beta_switch,beta_death=beta_death_high_effect,beta_cens=beta_cens),constant)
  scen_trt_sw_H1_low<-cbind(params_scenarios_grid(mu_W=mu_W,mu_L=mu_L,Sigma_W_L=Sigma_W_L,beta_prog=beta_prog,beta_switch=beta_switch,beta_death=beta_death_low_effect,beta_cens=beta_cens),constant)
  scen_trt_sw_H0_high<-cbind(params_scenarios_grid(mu_W=mu_W,mu_L=mu_L,Sigma_W_L=Sigma_W_L,beta_prog=beta_prog,beta_switch=beta_switch,beta_death=beta_death_H0_high_effect,beta_cens=beta_cens),constant)
  scen_trt_sw_H0_low<-cbind(params_scenarios_grid(mu_W=mu_W,mu_L=mu_L,Sigma_W_L=Sigma_W_L,beta_prog=beta_prog,beta_switch=beta_switch,beta_death=beta_death_H0_low_effect,beta_cens=beta_cens),constant)

  all_scen<-rbind(
    scen_trt_sw_H1_high,
    scen_trt_sw_H1_low,
    scen_trt_sw_H0_high,
    scen_trt_sw_H0_low
  )

  #additional specifics
  #if L has a pre-set effect on death, it must also have an effect on switching and progression:
  #if L has a pre-set effect on prog or on switch or on censoring, it must also have an effect on death:
  #if the mean trajectory of L is special (other than base pattern "1"), there should be full confounding
  for(i in 1:dim(all_scen)[1]) {
    case<-0

    if(all_scen$beta_death[[i]][5] != 0) case<-1
    if(all_scen$beta_prog[[i]][5] != 0 | all_scen$beta_switch[[i]][5] != 0 | all_scen$beta_cens[[i]][5] != 0 ) case<-2
    if(any(unlist(all_scen$mu_L[[i]])!=0)) case<-3

    if(case==1 | case==3) {
      all_scen$beta_prog[[i]][5]<-log(0.5)
      all_scen$beta_switch[[i]][5]<-log(1.5)
    }
    if(case==2 | case==3) all_scen$beta_death[[i]][5]<-log(0.5) #index 5 is effect of L

    #in th H0 scnenario with unequal trajectories of W, there should be no effect of W (index 3) on death
    if(!all(all_scen$mu_W[[i]][[1]]==all_scen$mu_W[[i]][[2]]) & all_scen$beta_death[[i]][6] == 0) all_scen$beta_death[[i]][3]<-0
    if(!all(all_scen$mu_L[[i]][[1]]==all_scen$mu_L[[i]][[2]]) & all_scen$beta_death[[i]][6] == 0) all_scen$beta_death[[i]][5]<-0

 }



  all_scen
}


#' Calculate true summary statistics for scenarios with delayed treatment effect
#'
#' @param Design Design data.frame for x
#' @param fixed_objects fixed objects not used for now
#'
#' @return For
#'
#' @export
#'
#' @describeIn generate_oncology  calculate true summary statistics for ...
#'
#' @examples
#' oncology_scenario_set_truevalues(oncology_scenario())
oncology_scenario_set_truevalues <- function(Design, fixed_objects = NULL) {
  oncology_scenario_set_truevalues_rowwise <- function(condition, cutoff_stats) {
    if (!"ev_soll" %in% colnames(condition)) {
      condition$ev_soll <- ceiling(((qnorm(1 - condition$alpha / 2) + qnorm(condition$power))
      / as.numeric(condition$beta_death[[1]][8]))^2 / condition$p_trt / (1 - condition$p_trt))
      condition
    } else {
      condition
    }
  }
  Design <- Design |>
    split(1:nrow(Design)) |>
    mapply(FUN = oncology_scenario_set_truevalues_rowwise, SIMPLIFY = FALSE)
  Design <- do.call(rbind, Design)
  Design
}
