library(tidyverse)
library(SimDesign)
library(survival)

source("generate_MACE.R")



#' Function that censor MACE events after discontinuation + buffer window
#' t_mace = T_disc + buffer is censored at one year
#'
#' @param data data in wide format, output of generate_mace function
#' @param buffer real non-negative value, duration of the buffer window
#'
#' @return data in wide format, with MACE events censored after the buffer window
#' @export
#'
#' @examples
#' toy_data <- generate_mace(Design[1,])
#' toy_data_trunc <- censor_buffer_window(toy_data,buffer)
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



#' Cox model without covariates
#' @param data data in wide format, output of censor_buffer window
#' 
#' @importFrom survival coxph
#' 
#' @return a coxph fit 
#' @export
#'
#' @examples
#' toy_data <- generate_mace(Design[1,])
#' toy_data_trunc <- censor_buffer_window(toy_data,buffer)
#' cox_model_nocov(toy_data_trunc)
cox_model_nocov <- function(data) {
  fit <- coxph(Surv(t_mace, event_mace) ~ A, data = data, id = ID)
  return(fit)
}

#' Cox model with covariates
#' @param data data in wide format, output of censor_buffer window
#'
#' @importFrom survival coxph
#' 
#' @return a coxph fit 
#' @export
#'
#' @examples
#' toy_data <- generate_mace(Design[1,])
#' toy_data_trunc <- censor_buffer_window(toy_data,buffer)
#' cox_model_cov(toy_data_trunc)
cox_model_cov <- function(data) {
  fit <- coxph(Surv(t_mace, event_mace) ~ A+X+Z, data = data, id = ID)
  return(fit)
}



#' Function that transforms data set with t_mace,event_mace, with one row per patients
#' into a data set with t_mace_start, t_mace_stop, event_mace, with discontinuation as
#' a time varying covariate
#'
#' @param data simulated data in wide format, output of generate_mace and censor_buffer_window
#'
#' @return data in long format for IPW analysis
#' @export
#'
#' @examples
#' toy_data <- generate_mace(Design[1,])
#' toy_data_trunc <- censor_buffer_window(toy_data,buffer)
#' toy_data_trunc_long <- data_process_start_stop(toy_data_trunc)
data_process_start_stop <- function(data) {
  disc_not_withdraw <- data$event_disc==1 & data$withdraw==0
  data$disc=0
  data2 <- data[disc_not_withdraw,]
  
  data$t_mace_start <- 0
  data$t_mace_stop <- ifelse(data$event_disc==1,data$t_disc,ifelse(data$event_mace==1,data$t_mace,1))
  
  data2$t_mace_start <- data2$t_disc
  data2$t_mace_stop <- data2$t_mace
  data2$disc <- 1
  
  return(rbind(data,data2)|>arrange(ID,t_mace_start))
}


#' Cox model with IPW weights with covariates
#'
#' @param data simulated data with long format, output of data_process_start_stop function
#' 
#' @importFrom survival coxph
#' @importFrom stats glm
#'
#' @return a coxph fit
#' @export
#'
#' @examples
#' toy_data <- generate_mace(Design[1,])
#' toy_data_trunc <- censor_buffer_window(toy_data,buffer)
#' toy_data_trunc_long <- data_process_start_stop(toy_data_trunc)
#' cox_model_ipw(toy_data_trunc_long)
cox_model_ipw <- function(data) {
  
  dat_withdraw <- data |> filter(event_disc==1,!duplicated(ID))
  fit_withdraw_0 <- glm(withdraw~X+Z,data=dat_withdraw|>filter(A==0),family = binomial(link = "logit"))
  fit_withdraw_1 <- glm(withdraw~X+Z,data=dat_withdraw|>filter(A==1),family = binomial(link = "logit"))
  
  
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




#' Cox model with IPW weights without covariates
#'
#' @param data simulated data with long format, output of data_process_start_stop function
#' 
#' @importFrom survival coxph
#' @importFrom stats glm
#' 
#' @return a coxph fit
#' @export
#'
#' @examples
#' toy_data <- generate_mace(Design[1,])
#' toy_data_trunc <- censor_buffer_window(toy_data,buffer)
#' toy_data_trunc_long <- data_process_start_stop(toy_data_trunc)
#' cox_model_ipw_nocov(toy_data_trunc_long)
cox_model_ipw_nocov <- function(data) {
  
  dat_withdraw <- data |> filter(event_disc==1,!duplicated(ID))
  fit_withdraw_0 <- glm(withdraw~X+Z,data=dat_withdraw|>filter(A==0),family = binomial(link = "logit"))
  fit_withdraw_1 <- glm(withdraw~X+Z,data=dat_withdraw|>filter(A==1),family = binomial(link = "logit"))
  
  
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


#' Function estimating all analysis functions on a simualted data set
#'
#' @param condition row of Design dataset
#' @param dat output from generate_mace function
#' @param fixed_objects list of parameters that are fixed across simulations
#'
#' @return named list with estimation metrics from the simulated replicate
#' @export
#'
#' @examples
#' for internal use by RunSimulation function
analyse_mace <- function(condition, dat, fixed_objects) {
  assumed_window <- condition$assumed_window
  
  trunc_data <- censor_buffer_window(dat,assumed_window)
  long_trunc_data <- data_process_start_stop(trunc_data)
  
  cox_nocov <- cox_model_nocov(trunc_data)
  cox_cov <- cox_model_cov(trunc_data)
  
  cox_nocov.trt <- coef(cox_nocov)["A"]
  cox_nocov.trt.se <- summary(cox_nocov)$coefficient[1,"se(coef)"]
  cox_nocov.HR <- exp(cox_nocov.trt)
  cox_nocov.HR.CI.l <- summary(cox_nocov)$conf.int[1,3]
  cox_nocov.HR.CI.u <- summary(cox_nocov)$conf.int[1,4]
  
  cox_cov.trt <- coef(cox_cov)["A"]
  cox_cov.trt.se <- summary(cox_cov)$coefficient[1,"se(coef)"]
  cox_cov.HR <- exp(cox_cov.trt)
  cox_cov.HR.CI.l <- summary(cox_cov)$conf.int[1,3]
  cox_cov.HR.CI.u <- summary(cox_cov)$conf.int[1,4]
  
  
  ipw_nocov <- cox_model_ipw_nocov(long_trunc_data)
  ipw_cov <- cox_model_ipw(long_trunc_data)
  
  ipw_nocov.trt <- coef(ipw_nocov)["A"]
  ipw_nocov.trt.se <- summary(ipw_nocov)$coefficient[1,"robust se"]
  ipw_nocov.HR <- exp(ipw_nocov.trt)
  ipw_nocov.HR.CI.l <- summary(ipw_nocov)$conf.int[1,3]
  ipw_nocov.HR.CI.u <- summary(ipw_nocov)$conf.int[1,4]
  
  ipw_cov.trt <- coef(ipw_cov)["A"]
  ipw_cov.trt.se <- summary(ipw_cov)$coefficient[1,"robust se"]
  ipw_cov.HR <- exp(ipw_cov.trt)
  ipw_cov.HR.CI.l <- summary(ipw_cov)$conf.int[1,3]
  ipw_cov.HR.CI.u <- summary(ipw_cov)$conf.int[1,4]
  
  res <- nc(cox_nocov.trt=cox_nocov.trt,
            cox_nocov.trt.se=cox_nocov.trt.se,
            cox_nocov.HR=cox_nocov.HR,
            cox_nocov.HR.CI.l=cox_nocov.HR.CI.l,
            cox_nocov.HR.CI.u=cox_nocov.HR.CI.u,
            
            cox_cov.trt=cox_cov.trt,
            cox_cov.trt.se=cox_cov.trt.se,
            cox_cov.HR=cox_cov.HR,
            cox_cov.HR.CI.l=cox_cov.HR.CI.l,
            cox_cov.HR.CI.u=cox_cov.HR.CI.u,
            
            ipw_nocov.trt=ipw_nocov.trt,
            ipw_nocov.trt.se=ipw_nocov.trt.se,
            ipw_nocov.HR=ipw_nocov.HR,
            ipw_nocov.HR.CI.l=ipw_nocov.HR.CI.l,
            ipw_nocov.HR.CI.u=ipw_nocov.HR.CI.u,
            
            ipw_cov.trt=ipw_cov.trt,
            ipw_cov.trt.se=ipw_cov.trt.se,
            ipw_cov.HR=ipw_cov.HR,
            ipw_cov.HR.CI.l=ipw_cov.HR.CI.l,
            ipw_cov.HR.CI.u=ipw_cov.HR.CI.u
            )
  
  res
  
}


#' Function that calculates all simulation metrics for each scenario 
#'
#' @param condition row of Design dataset
#' @param results output from analyse_mace function
#' @param fixed_objects list of parameters that are fixed across simulations

#'
#' @return named list of summarized results from the given Design row
#' @export
#'
#' @examples
#' internal use by RunSimulation
summarise_mace <- function(condition, results, fixed_objects) {
  true_trt <- as.numeric(condition$true_trt)
  
  with(results, {
    
list_results <- c(    
    bias_cox_nocov_trt = mean(cox_nocov.trt-true_trt),
    bias_cox_cov_trt = mean(cox_cov.trt-true_trt),
    bias_ipw_nocov_trt = mean(ipw_nocov.trt-true_trt),
    bias_ipw_cov_trt = mean(ipw_cov.trt-true_trt),
    
    mean_se_trt_cox_nocov = mean(cox_nocov.trt.se),
    mean_se_trt_cox_cov = mean(cox_cov.trt.se),
    mean_se_trt_ipw_nocov = mean(ipw_nocov.trt.se),
    mean_se_trt_ipw_cov = mean(ipw_cov.trt.se),
    
    emp_se_trt_cox_nocov = sd(cox_nocov.trt),
    emp_se_trt_cox_cov = sd(cox_cov.trt),
    emp_se_trt_ipw_nocov = sd(ipw_nocov.trt),
    emp_se_trt_ipw_cov = sd(ipw_cov.trt),
    
    mean_HR_cox_nocov = mean(cox_nocov.HR),
    mean_HR_cox_cov = mean(cox_cov.HR),
    mean_HR_ipw_nocov = mean(ipw_nocov.HR),
    mean_HR_ipw_cov = mean(ipw_cov.HR),
    
    mean_HR_CI_width_cox_nocov = mean(cox_nocov.HR.CI.u-cox_nocov.HR.CI.l),
    mean_HR_CI_width_cox_cov = mean(cox_cov.HR.CI.u-cox_cov.HR.CI.l),
    mean_HR_CI_width_ipw_nocov = mean(ipw_nocov.HR.CI.u-ipw_nocov.HR.CI.l),
    mean_HR_CI_width_ipw_cov = mean(ipw_cov.HR.CI.u-ipw_cov.HR.CI.l),
    
    power_cox_nocov = mean((cox_nocov.trt^2)/(cox_nocov.trt.se^2)>qchisq(0.95,1)),
    power_cox_cov = mean((cox_cov.trt^2)/(cox_cov.trt.se^2)>qchisq(0.95,1)),
    power_ipw_nocov = mean((ipw_nocov.trt^2)/(ipw_cov.trt.se^2)>qchisq(0.95,1)),
    power_ipw_cov = mean((ipw_cov.trt^2)/(ipw_cov.trt.se^2)>qchisq(0.95,1))
  
  )
return(list_results)
  })
  
}

Design <- assumptions_mace()

complete_design <- true_trt_mace(Design)

results <- runSimulation(
  complete_design,
  replications = 1000,
  generate = generate_mace,
  analyse = analyse_mace,
  summarise = summarise_mace,
  fixed_objects = c()
)

### example of results for power and bias

results[,grep("power",colnames(results))]
results[,grep("bias",colnames(results))]
