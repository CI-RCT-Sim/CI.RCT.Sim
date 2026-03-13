library(tidyverse)
library(SimDesign)
library(survival)
library(survminer)

source("generate_MACE.R")



###     censor MACE events after discontinuation + buffer window
###     t_mace = T_disc + buffer is censored at one year

censor_buffer_window <- function(data,buffer) {
  disc_not_withdraw <- data$event_disc==1 & data$withdraw==0
  
  data[disc_not_withdraw,"event_mace"] <- ifelse(data[disc_not_withdraw,"t_mace"]-data[disc_not_withdraw,"t_disc"]<buffer,
                                                 data[disc_not_withdraw,"event_mace"],
                                                 0)
  
  data[disc_not_withdraw,"t_mace"] <- ifelse(data[disc_not_withdraw,"t_mace"]-data[disc_not_withdraw,"t_disc"]<buffer,
                                                 data[disc_not_withdraw,"t_mace"],
                                                 pmin(1,data[disc_not_withdraw,"t_disc"]+buffer))
  return(data)
}




cox_model_nocov <- function(data) {
  fit <- coxph(Surv(t_mace, event_mace) ~ A, data = data, id = ID)
  return(fit)
}

cox_model_cov <- function(data) {
  fit <- coxph(Surv(t_mace, event_mace) ~ A+X+Z, data = data, id = ID)
  return(fit)
}


### function that transforms data set with t_mace,event_mace, with one row per patients
### into a data set with t_mace_start, t_mace_stop, event_mace, with discontinuation as
### a time varying covariate

data_process_start_stop <- function(data) {
  disc_not_withdraw <- data$event_disc==1 & data$withdraw==0
  data$disc=0
  data2 <- data[disc_not_withdraw,]
  
  data$tstart <- 0
  data$tstop <- ifelse(data$event_disc==1,data$t_disc,ifelse(data$event_mace==1,data$t_mace,1))
  
  data2$tstart <- data2$t_disc
  data2$tstop <- data2$t_mace
  data2$disc <- 1
  
  return(rbind(data,data2)%>%arrange(ID,tstart))
}


cox_model_ipw <- function(data,buffer) {
  
  dat_withdraw <- data %>% filter(event_disc==1,!duplicated(ID))
  fit_withdraw_0 <- glm(withdraw~X+Z,data=dat_withdraw%>%filter(A==0),family = binomial(link = "logit"))
  fit_withdraw_1 <- glm(withdraw~X+Z,data=dat_withdraw%>%filter(A==1),family = binomial(link = "logit"))
  
  
  p0 <- predict(fit_withdraw_0, newdata = data, type = "response")
  p1 <- predict(fit_withdraw_1, newdata = data, type = "response")
  lp <- c()
  for (i in 1:dim(data)[1]) {
    if (data$disc[i]==1) {
      if (data[i,"A"]==0) {
        lp <- c(lp,1-p0[i])
      }else {
        lp <- c(lp,1-p1[i])
      }
      
    } else {
      lp <- c(lp,1)
    } 
  } 
  data$IPW <- 1/lp
  (coxph(Surv(t_mace_start,t_mace_stop,event_mace)~A+X+Z,data=data,id=ID,weights=IPW))
  
}

cox_model_ipw_nocov <- function(data,buffer) {
  
  dat_withdraw <- data %>% filter(event_disc==1,!duplicated(ID))
  fit_withdraw_0 <- glm(withdraw~X+Z,data=dat_withdraw%>%filter(A==0),family = binomial(link = "logit"))
  fit_withdraw_1 <- glm(withdraw~X+Z,data=dat_withdraw%>%filter(A==1),family = binomial(link = "logit"))
  
  
  p0 <- predict(fit_withdraw_0, newdata = data, type = "response")
  p1 <- predict(fit_withdraw_1, newdata = data, type = "response")
  lp <- c()
  for (i in 1:dim(data)[1]) {
    if (data$disc[i]==1) {
      if (data[i,"A"]==0) {
        lp <- c(lp,1-p0[i])
      }else {
        lp <- c(lp,1-p1[i])
      }
      
    } else {
      lp <- c(lp,1)
    } 
  } 
  data$IPW <- 1/lp
  (coxph(Surv(t_mace_start,t_mace_stop,event_mace)~A,data=data,id=ID,weights=IPW))
  
}


#### Toy example ####
assumed_window <- as.numeric(Design[1,"assumed_window"])

toy_data <- Generate(Design[1,])
toy_data_buffer <- censor_buffer_window(toy_data,assumed_window)

cox_model_cov(toy_data_buffer)

toy_data_start <- data_process_start_stop(toy_data_buffer)
cox_model_ipw(toy_data_start,buffer=assumed_window)
