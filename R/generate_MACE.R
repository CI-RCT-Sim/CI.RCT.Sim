library(tidyverse)
library(SimDesign)


baseline <- list(

  ########## MACE parameters ##########
  beta_mace_0 = c(log(-log(1-0.1))),
  beta_mace_trt_before = c(log(1.5)),
  beta_mace_trt_buffer = c(log(1.5)),
  beta_mace_post = c(0),
  beta_mace_X = c(log(2)),
  beta_mace_Z = c(log(2)),
  beta_mace_L = c(0),

  ########## Discontinuation parameters ##########
  beta_disc_0 = c(log(-log(1-0.5))),
  beta_disc_trt = c(log(2)),
  beta_disc_X = c(log(2)),
  beta_disc_Z = c(log(2)),
  beta_disc_L = c(0),

  ########## Withdrawal parameters ##########
  beta_wd_0 = c(log(0.5/(1-0.5))),
  beta_wd_trt = c(log(2)),
  beta_wd_X = c(log(2)),
  beta_wd_Z = c(log(2)),
  beta_wd_L = c(0),

  ########## Buffer window ##########
  true_window = c(1/12),
  assumed_window = c(1/12)
)


### mace parameters  ###
further1 <- baseline%>%replace("beta_mace_0",log(-log(1-0.3))) # higher beta_mace_0
further2 <- baseline%>%replace("beta_mace_trt_before",0) # no TE
further3 <- baseline%>%replace("beta_mace_trt_before",log(2)) # greater TE
further4 <- baseline%>%replace("beta_mace_trt_buffer",0) # no effect in buffer
further5 <- baseline%>%replace("beta_mace_trt_buffer",log(1.25)) # lower effect in buffer
further6 <- baseline%>%replace(c("beta_mace_X","beta_mace_Z"),c(0,0)) # no covariate effects
further7 <- baseline%>%replace(c("beta_mace_L"),log(2)) # confounder covariate effect

### disc parameters ###
further8 <- baseline%>%replace("beta_disc_0",log(-log(1-0.2))) # lower beta_disc_0
further9 <- baseline%>%replace(c("beta_disc_X","beta_disc_Z"),c(0,0)) # no covariate effects
further10 <- baseline%>%replace(c("beta_disc_L"),log(2)) # confounder covariate effect

## withdrawal parameters ###
further11 <- baseline%>%replace("beta_wd_0",log(0.75/(1-0.75))) # greater withdrawal
further12 <- baseline%>%replace(c("beta_wd_X","beta_wd_Z"),c(0,0)) # no covariate effects
further13 <- baseline%>%replace(c("beta_wd_L"),log(2)) # confounder covariate effect

### buffer window ###
further14 <- baseline%>%replace(c("assumed_window"),0) # lower assumed window
further15 <- baseline%>%replace(c("assumed_window"),2/12) # greater assumed window


OFATT <- bind_rows(baseline,further1,further2,further3,further4,further5,
                   further6,further7,further8,further9,further10,further11,
                   further12,further13,further14,further15
                   )


OFATT$logHRassumed <- OFATT$beta_mace_trt_before              ### log HR assumed fixed to beta_mace_trt_before for H1

OFATT[OFATT$beta_mace_trt_before==0,"logHRassumed"] <- log(1.5) ### for H0, fixed to either log(1.5) or log(2)

OFATT2 <- OFATT%>%filter(beta_mace_trt_before==0)%>%
  mutate(logHRassumed=log(2))
design_list <- as.list(rbind(OFATT,OFATT2))
design_list$fully.crossed <- F

Design <- do.call(createDesign,design_list)

############# function to generate one TTE MACE #######################

###     vectorised function takes the data from one patient, and returns
###     t_mace, event_mace, t_disc censored at one year. At this stage we don't
###     censor MACE events outside the buffer window
###     with no_withdraw=T, there is no withdrawal for the estimation of the true TE

# cov_df <- data.frame(X=X,Z=Z,L=L,A=A,hazard_disc=hazard_disc,logit_withdraw=logit_withdraw,beta_tot_before=beta_tot_before,beta_tot_buffer=beta_tot_buffer,beta_tot_after=beta_tot_after,ID=ID)

generate_all_tte2 <- function(cov_df, condition=condition, no_withdraw=F) {
  with(cov_df, {
    N <- length(X)
    true_window <- condition$true_window

    ## checking if we should consider withdrawal
    if (no_withdraw==T) {
      prob_withdraw <- 0
    } else {
      prob_withdraw <- exp(logit_withdraw)/(1+exp(logit_withdraw))

    }
    #################### Simulate discontinuation times ####################
    u_disc <- runif(N,0,1)
    T_disc <- -log(u_disc)/cov_dv$hazard_disc

    t_0 <- T_disc
    t_1 <- T_disc+true_window

    ### Simulate MACE times, with a time-varying treatment effect with two breakpoints
    ### t_0 = T_disc and T_1 = T_disc + buffer

    u_mace <- runif(N,0,1)
    m_logu <- -log(u_mace)
    if ((m_logu <= t_0*exp(beta_tot_before))) {
      T_mace <- m_logu/(exp(beta_tot_before))
    } else {
      if ((m_logu > t_0*exp(beta_tot_before) & m_logu<= t_0*exp(beta_tot_before) +(t_1-t_0)*exp(beta_tot_buffer))) {
        T_mace <-(m_logu-t_0*exp(beta_tot_before)+(t_0)*exp(beta_tot_buffer))/exp(beta_tot_buffer)
      } else {
        T_mace <- t_1+(m_logu -t_0*exp(beta_tot_before)-(t_1-t_0)*exp(beta_tot_buffer))/(exp(beta_tot_after))
      }
    }


    data <- data.frame(ID=ID,X=X,Z=Z,A=A,T_disc=T_disc,T_mace=T_mace,u_mace=u_mace,prob_withdraw,m_logu=m_logu)
    data$disc_before_mace <- T_disc<T_mace
    data$censor_mace = T_mace>1
    data$censor_disc = T_disc>1

    dat_i <- data
    dat_temp <- dat_i%>%select(-T_mace)     # MC comment 3
    data_formatted <- c()

    ### simulate withdrawal
    u_withdraw <- runif(N,0,1)
    is_withdraw <- u_withdraw<prob_withdraw

    if (is_withdraw & (dat_i$T_mace>dat_i$T_disc)) {
      ### censoring only if disc is before MACE and patient withdrew, with
      ### administrative censoring at one year
      dat_temp$t_mace <- pmin(dat_i$T_disc,1)
      dat_temp$event_mace <- c(0)
      dat_temp$censor  <- 1
      dat_temp$withdraw  <- dat_i$T_disc<1
      dat_temp$t_disc <- pmin(dat_i$T_disc,1)
      dat_temp$event_disc <- dat_i$T_disc<1

    } else {
      ### otherwise, record MACE event, with administrative censoring at one year
      dat_temp$t_mace <- pmin(dat_i$T_mace,1)
      dat_temp$event_mace <- ifelse(dat_i$T_mace>1,0,1)
      dat_temp$censor  <- ifelse(dat_i$T_mace>1,1,0)
      dat_temp$withdraw  <- 0
      dat_temp$t_disc <- ifelse(dat_i$T_mace<dat_i$T_disc,pmin(dat_i$T_mace,1),pmin(dat_i$T_disc,1))
      dat_temp$event_disc <- ifelse(dat_i$T_mace<dat_i$T_disc,0,dat_i$T_disc<1)
    }
    data_formatted <- dat_temp
    data_formatted$T_mace <- T_mace
    return(data_formatted)
  })
}


Generate <- function(condition, fixed_objects) {

  ########## MACE parameters ##########
  beta_mace_0 = condition$beta_mace_0
  beta_mace_trt_before = condition$beta_mace_trt_before
  beta_mace_trt_buffer = condition$beta_mace_trt_buffer
  beta_mace_post = condition$beta_mace_post
  beta_mace_X = condition$beta_mace_X
  beta_mace_Z = condition$beta_mace_Z
  beta_mace_L = condition$beta_mace_L

  ########## Discontinuation parameters ##########
  beta_disc_0 = condition$beta_disc_0
  beta_disc_trt = condition$beta_disc_trt
  beta_disc_X = condition$beta_disc_X
  beta_disc_Z = condition$beta_disc_Z
  beta_disc_L = condition$beta_disc_L

  ########## Withdrawal parameters ##########
  beta_wd_0 = condition$beta_wd_0
  beta_wd_trt = condition$beta_wd_trt
  beta_wd_X = condition$beta_wd_X
  beta_wd_Z = condition$beta_wd_Z
  beta_wd_L = condition$beta_wd_L

  ########## Buffer window ##########
  true_window = condition$true_window
  assumed_window = condition$assumed_window

  ########## Assumed HR #############
  logHR_assumed = condition$logHRassumed

  #################### Compute sample size  ####################
  e_mace = 4*(qnorm(1-0.05/2)+qnorm(1-0.2))^2/(logHR_assumed^2)  ### Schoenfeld’s formula
  # MC: should this be e_mace <- ((qnorm(1-0.05/2) + qnorm(1-0.2))/logHR_assumed)^2 * 4   ?
  sample_size = floor(e_mace/(1-exp(-exp(beta_mace_0))))  ####### sample size calculation

  #################### Simulate covariates  ####################
  X <- rnorm(sample_size,0,1)
  Z <- rnorm(sample_size,0,1)
  L <- rnorm(sample_size,0,1)
  A <- sample(0:1,size=sample_size,replace=T,prob=c(0.5,0.5))

  hazard_disc <- exp(beta_disc_0+(0.5-A)*beta_disc_trt+X*beta_disc_X+Z*beta_disc_Z+L*beta_disc_L)
  logit_withdraw <- beta_wd_0+X*beta_wd_X+Z*beta_wd_Z+L*beta_wd_L+beta_wd_trt*(1-A)
  beta_fix_mace <- beta_mace_0+X*beta_mace_X+Z*beta_mace_Z+L*beta_mace_L

  ## Defining time-varying treatment effect after discontinuation
  beta_tot_before <- beta_fix_mace+(1-A)*beta_mace_trt_before
  beta_tot_buffer <- beta_fix_mace+(1-A)*beta_mace_trt_buffer
  beta_tot_after <- beta_fix_mace
  cov_df <- data.frame(X=X,Z=Z,L=L,A=A,hazard_disc=hazard_disc,logit_withdraw=logit_withdraw,beta_tot_before=beta_tot_before,beta_tot_buffer=beta_tot_buffer,beta_tot_after=beta_tot_after,ID=1:sample_size)
  #################### Simulate TTE  ####################

  res <- generate_all_tte2(cov_df,condition=condition,no_withdraw=FALSE)

  dat_f <- res %>%
    select(-prob_withdraw,-u_mace,-m_logu,-censor_mace,-censor_disc,-disc_before_mace,-prob_withdraw)

  return(dat_f)
}

calc_true_trt <- function(condition,fixed_objects) {

  ########## MACE parameters ##########
  beta_mace_0 = condition$beta_mace_0
  beta_mace_trt_before = condition$beta_mace_trt_before
  beta_mace_trt_buffer = condition$beta_mace_trt_buffer
  beta_mace_post = condition$beta_mace_post

  ########## Discontinuation parameters ##########
  beta_disc_0 = condition$beta_disc_0
  beta_disc_trt = condition$beta_disc_trt

  ########## Buffer window ##########
  true_window = condition$true_window
  assumed_window = condition$assumed_window

  ########## Assumed HR #############
  logHR_assumed = condition$logHRassumed

  count_events_mace <- 0
  X <- Z <- L <- rep(0,1000000)
  A <-sample(0:1,1000000,prob = c(0.5,0.5),replace=T)
  ID <- 1:1000000
  # cov_df <- data.frame(X=rep(0,1000000),Z=rep(0,1000000),L=rep(0,1000000),A=sample(0:1,1000000,prob = c(0.5,0.5),replace=T),)

  hazard_disc <- exp(beta_disc_0+(0.5-A)*beta_disc_trt)
  logit_withdraw <- 0
  beta_fix_mace <- beta_mace_0

  ## Defining time-varying treatment effect after discontinuation
  beta_tot_before <- beta_fix_mace+(1-A)*beta_mace_trt_before
  beta_tot_buffer <- beta_fix_mace+(1-A)*beta_mace_trt_buffer
  beta_tot_after <- beta_fix_mace

  cov_df <- data.frame(X=X,Z=Z,L=L,A=A,hazard_disc=hazard_disc,logit_withdraw=logit_withdraw,beta_tot_before=beta_tot_before,beta_tot_buffer=beta_tot_buffer,beta_tot_after=beta_tot_after,ID=ID)

  dat_bind_bl <- generate_all_tte2(cov_df, condition=condition, no_withdraw=T)
  count_events_mace <- sum(dat_bind_bl$event_mace)
  dat_bind <- c()

  dat <- rbind(dat_bind_bl,bind_rows(dat_bind))
  dat <- transform_disc_tte(dat)%>%select(-censor)
  dat2 <- censor_buffer_window(dat,true_window)
  fit <- coxph(Surv(t_mace, event_mace) ~ A, data = dat2, id = ID)
  true_trt <- c()
  true_trt[[1]] <- fit
  true_trt[[2]] <- dat2
  return(true_trt)
}
