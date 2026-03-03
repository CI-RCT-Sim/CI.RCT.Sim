library(tidyverse)
library(SimDesign)
library(survival)
library(survminer)

source("data_generation_scenario_4.R")



###     censor MACE events after discontinuation + buffer window
###     t_mace = T_disc + buffer can be greater than one year while T_disc<1

within_buffer_apply <- function(ID,X,Z,A,t_disc,t_mace,event_mace,T_mace,event_disc,withdraw,buffer,...) {
  if (t_disc< t_mace) { ### if disc before mace
    if (between(T_mace,t_disc,t_disc+buffer)) { ### if  within buffer
      t_mace <-T_mace 
      event_mace<- 1
    } else { ### not within buffer
      t_mace <- (t_disc+buffer)
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


cox_model_ipw <- function(data,buffer) {
  
  fit_disc_0 <- coxph(Surv(t_disc, event_disc) ~ X+Z, data = data%>%filter(A==0), id = ID)
  fit_disc_1 <- coxph(Surv(t_disc, event_disc) ~ X+Z, data = data%>%filter(A==1), id = ID)
  
  sf_0 <- survfit(fit_disc_0,newdata = data)
  sf_1 <- survfit(fit_disc_1,newdata = data)
  
  data$surv_prob <- sapply(1:nrow(data), function(i) {
    if (data$A[i]==0) {
      return(summary(sf_0[i], times = data$t_mace[i], extend = TRUE)$surv)
    } else {
      return(summary(sf_1[i], times = data$t_mace[i], extend = TRUE)$surv)
    }
    
  })
  
  data$surv_prob2 <- sapply(1:nrow(data), function(i) {
    if (data$A[i]==0) {
      return(summary(sf_0[i], times = max(data$t_mace[i]-buffer,0), extend = TRUE)$surv)
    } else {
      return(summary(sf_1[i], times = max(data$t_mace[i]-buffer,0), extend = TRUE)$surv)
    }
    
  })
  dat_withdraw <- data %>% filter(event_disc==1)
  
  fit_withdraw_0 <- glm(withdraw~X+Z,data=dat_withdraw%>%filter(A==0),family = binomial(link = "logit"))
  fit_withdraw_1 <- glm(withdraw~X+Z,data=dat_withdraw%>%filter(A==1),family = binomial(link = "logit"))
  
  
  p0 <- predict(fit_withdraw_0, newdata = data, type = "response")
  p1 <- predict(fit_withdraw_1, newdata = data, type = "response")
  
  lp <- c()
  for (i in 1:dim(data)[1]) {
    if (data[i,"A"]==0) {
      lp <- c(lp,p0[i])
    }else {
      lp <- c(lp,p1[i])
    }
    
  } 
  data$prob_withdraw <- lp
  
  eps <- 0.00001
  data$prob_withdraw <- pmin(pmax(data$prob_withdraw, eps), 1 - eps) 
  p_miss2 <- (data$surv_prob2-data$surv_prob)*(data$prob_withdraw)
  p_miss2 <- pmin(pmax(p_miss2, eps), 1 - eps) 
  
  data$IPW <- 1 /(1-p_miss2)
  
  fit_mace_nocov <- coxph(Surv(t_mace, event_mace) ~ A,data=data,weights=IPW,id=ID)
  fit_mace_cov <- coxph(Surv(t_mace, event_mace) ~ A+X+Z,data=data,weights=IPW,id=ID)  
  fit_final <- c()
  fit_final[[1]] <- fit_mace_nocov
  fit_final[[2]] <- fit_mace_cov
  
  return(fit_final)
}

#### Toy example ####
assumed_window <- as.numeric(Design[1,"assumed_window"])

toy_data <- Generate(Design[1,])
toy_data_buffer <- censor_buffer_window(toy_data,assumed_window)

cox_model_cov(toy_data_buffer)
cox_model_ipw(toy_data_buffer,buffer=assumed_window)
