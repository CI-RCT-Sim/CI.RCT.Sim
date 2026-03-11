library(tidyverse)
library(SimDesign)
library(survival)
library(survminer)

source("data_generation_scenario_4.R")



###     censor MACE events after discontinuation + buffer window
###     t_mace = T_disc + buffer is censored at one year

within_buffer_apply <- function(ID,X,Z,A,t_disc,t_mace,event_mace,T_mace,event_disc,withdraw,buffer,...) {
  if (t_disc< t_mace) { ### if disc before mace
    if (between(t_mace,t_disc,t_disc+buffer)) { ### if  within buffer
      t_mace <-t_mace 
      event_mace <- event_mace
    } else { ### not within buffer
      t_mace <- min(t_disc+buffer,1)
      event_mace <- 0
    }
  } 
  return(data.frame(ID=ID,X=X,Z=Z,A=A,t_disc=t_disc,t_mace=t_mace,event_mace=event_mace,event_disc=event_disc,withdraw=withdraw))
}

censor_buffer_window <- function(data,buffer) {
  return(pmap(data,within_buffer_apply,buffer=buffer)%>%bind_rows)
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
  df_final <- c()
  for (i in 1:dim(data)[1]) {
    dat_i <- data[i,]
    dat_i2 <- c() 
    dat_i$disc <- 0
    if (dat_i$event_disc==1) {
      dat_i$t_mace_start <- 0
      dat_i$t_mace_stop <- dat_i$t_disc
      dat_i$event_mace <- 0
      if (dat_i$withdraw==0) {
        dat_i2 <- data[i,]
        dat_i2$t_mace_start <- dat_i$t_disc
        dat_i2$t_mace_stop <- dat_i2$t_mace
        dat_i2$event_mace <- dat_i2$event_mace
        dat_i2$disc <- 1
      }
    } else {
      dat_i$t_mace_start <- 0
      dat_i$t_mace_stop <- dat_i$t_mace
      
    }
    df_final[[i]] <- rbind(dat_i,dat_i2)
  }
  return(bind_rows(df_final))
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
