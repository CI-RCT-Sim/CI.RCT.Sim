#Test

rm(list=ls())
#setwd("/home/robin/EMA_causal_16May_node2/CI.RCT.Sim")

devtools::load_all()
library(CI.RCT.Sim)
library(parallel)
library(survival)

sim_parameters <- oncology_scenario() |>
  oncology_scenario_set_truevalues()
#is(sim_parameters)
seed<-round(runif(1)*1000) #950, 400, 51, 871, 459
set.seed(seed)
A<-sim_parameters[1,] #2 oder 28
condition<-A
data<-generate_oncology(A,fixed_objects = list(allow_switch = TRUE, logHR_assumed = NULL, allow_random_cens =TRUE))


sum(data$ev)
sum(data$prog_ev)
sum(data$switch)
source("scripts/deskriptiv1-22utf.R")
head(data)
#save(data,file="results5/example_data_scen_2.RData")
destab<-make.table(varlist=c("X_0","W_0","X_2BL","W_2BL","ev","prog_ev","random_cens","switch","event_time"),gruppe="trt",r=2,r.p=4,r.proz=2,z.grenz=10,tab.data=data,aov.enforce=FALSE,my.test=NULL,testname="my",enforce.test=FALSE,nonpar=FALSE,collapse=FALSE,chisq.enforce=FALSE,fisher.enforce=FALSE,marg=1,notest=FALSE,rm_NA_only=TRUE)
destab
destab_median<-make.table(varlist=c("X_0","W_0","X_2BL","W_2BL","ev","prog_ev","random_cens","switch","event_time"),gruppe="trt",r=2,r.p=4,r.proz=2,z.grenz=10,tab.data=data,aov.enforce=FALSE,my.test=NULL,testname="my",enforce.test=FALSE,nonpar="range",collapse=FALSE,chisq.enforce=FALSE,fisher.enforce=FALSE,marg=1,notest=FALSE,rm_NA_only=TRUE)
destab_median
data$calendar_end_of_study[1]

library(survminer)


#summary(mod_ipw)
#mod<-mod_ipw
modsum<-function(mod) {
  HR<-exp(coef(mod)[1])
  CI <- exp(confint(mod)[1, ])
  smr <- summary(mod)
  p <- smr$coef[1, "Pr(>|z|)"]
  p<-ifelse(HR<=1,p/2,1-p/2)
  if(any(grepl("robust se",colnames(smr$coef)))) {
    SE <- smr$coef[1, "robust se"]
  } else {
    SE <- smr$coef[1, "se(coef)"]
  }
  data.frame(HR = HR, SElogHR = SE, low = CI[[1]],up = CI[[2]],  p = p)
}

#Treatment policy
mod_itt <- coxph(Surv(time = event_time, event = ev) ~ trt + X_0 + W_0, data = data)
itt_res<-modsum(mod_itt)


#Censored
data_cens<-data
set_cens <- data$trt == 0 & data$switch == 1
data_cens$event_time[set_cens] <- data$prog_time[set_cens]
data_cens$ev[set_cens] <- 0
mod_cens <- coxph(Surv(time = event_time, event = ev) ~ trt + X_0 + W_0, data = data_cens)
cens_res<-modsum(mod_cens)

#IPW
set<-data$trt == 0 & data$prog_ev == 1
wmod <- glm(switch ~ X_2BL + W_2BL, family = binomial, data = data, subset = set)
pred <- predict(wmod, type = "response")
pred_a <- ifelse(data$switch[set] == 1, pred, 1 - pred)
data$w <- 1
data$w[set] <- 1 / pred_a

sdat <- tmerge(data1 = data, data2 = data, id = id, tstop = event_time)
sdat <- tmerge(data1 = sdat, data2 = data, id = id, death_event = event(event_time, ev))
sdat <- tmerge(data1 = sdat, data2 = data, id = id, PD = tdc(prog_time))
sdat$w[sdat$PD == 0 | sdat$trt == 1] <- 1
#sdat$w<-0.5*sdat$w
#sdat$w<-dim(data)[1]*sdat$w/sum(sdat$w)
#ptrt<-mean(data$trt)
#sdat$w[sdat$trt==1]<-ptrt*sdat$w[sdat$trt==1]
#sdat$w[sdat$trt==0]<-(1-ptrt)*sdat$w[sdat$trt==0]

remo<-sdat$trt == 0 & sdat$PD == 1 & sdat$switch == 1
dat_a<-sdat[!remo,]
mod_ipw<-coxph(Surv(time=tstart,time2=tstop,event=death_event)~trt+X_0+W_0,data=dat_a,weights=w,robust=TRUE,id=id)
ipw_res<-modsum(mod_ipw)

boxplot(dat_a$w[dat_a$PD==1 & dat_a$trt==0])
summary(dat_a$w[dat_a$PD==1 & dat_a$trt==0])

#dat_a[dat_a$id==19,]
#dat_a$id[dat_a$trt==0]
#RPSFTM
prep_data_RPSFTM_fun <- function(data) {
  data$time_on_trt <- 0
  data$time_on_trt[data$trt == 1] <- data$event_time[data$trt == 1]
  set_sw <- data$trt == 0 & data$switch == 1
  data$time_on_trt[set_sw] <- data$event_time[set_sw] - data$prog_time[set_sw]
  data$time_on_trt_relative <- data$time_on_trt / data$event_time
  data$max_FU <- data$calendar_end_of_study - data$calendar_start_time
  # for admin censoring, max_FU is the event time. The numeric calclation may cause minimal differences
  # so
  set_admin_cens <- data$max_FU < data$event_time
  data$max_FU[set_admin_cens] <- data$event_time[set_admin_cens]
  data
}
data_rpsftm <- prep_data_RPSFTM_fun(data)
head(data_rpsftm)
data_rpsftm$max_FU

RPS <- rpsftm(
  data = data_rpsftm, id = "id", time = "event_time", event = "ev", treat = "trt",
  base_cov = c("X_0", "W_0"),
  rx = "time_on_trt_relative",
  psi_test = "phreg",
  alpha = 0.05,
  censor_time = "max_FU",
  autoswitch = TRUE,
  recensor = FALSE,
  gridsearch = FALSE,
  root_finding = "bisection",
  boot = TRUE,
  n_boot = 1000,
  nthreads = 1,
  seed = 851
)
RPS
plot(RPS)
RPS_res<-data.frame(HR = RPS$hr,SElogHR = sd(log(RPS$hr_boots),na.rm=TRUE),low = RPS$hr_CI[1],up = RPS$hr_CI[2],p = RPS$pvalue)
head(RPS$data_outcome)



#km_TSE<-survfit(Surv(time=t_star,event=d_star)~trt,data=TSE$data_outcome)
#plot(km_TSE)


#TSE
TSE <- tsesimp(
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
  alpha = 0.05,
  recensor = FALSE,
  swtrt_control_only = TRUE,
  offset = 0,
  boot = TRUE,
  n_boot = 1000,
  nthreads = 1,
  seed = 761
)

TSE_res<-data.frame(HR = TSE$hr,SElogHR = sd(log(TSE$hr_boots),na.rm=TRUE),low = TSE$hr_CI[1],up = TSE$hr_CI[2],p = TSE$pvalue)

###


#gformula

gform<-analyse_oncology_gformula_unequal(B=20, reps=10,n_ev_cutoff_no_bootstrap=1000,return_data=TRUE)(A,data)
gform_res<-as.data.frame(gform[1:5])



km_OS<-survfit(Surv(time=event_time,event=ev)~trt,data=data)
time_unit<-"Years"
OS_plot_rt<-ggsurvplot(km_OS,data=data,ylab="Overall survival",xlab=time_unit,risk.table=TRUE)
OS_plot_rt

km_OS_2BL<-survfit(Surv(time=event_time-prog_time,event=ev)~switch,data=data[data$trt==0 & data$prog_ev==1,])
OS_plot_2BL<-ggsurvplot(km_OS_2BL,data=data[data$trt==0 & data$prog_ev==1,],ylab="Overall survival",xlab=time_unit,risk.table=TRUE)
OS_plot_2BL

OS_plot<-ggsurvplot(km_OS,data=data,ylab="Overall survival",xlab=time_unit,title="ITT")
km_ipw<-survfit(Surv(time=tstart,time2=tstop,event=death_event)~trt,data=dat_a,weights=w)
OS_plot_ipw<-ggsurvplot(km_ipw,data=dat_a,ylab="Overall survival",xlab=time_unit,title="IPW")
km_rps<-survfit(Surv(time=t_star,event=d_star)~trt,data=RPS$data_outcome)
OS_plot_rps<-ggsurvplot(km_rps,ylab="Overall survival",xlab=time_unit,title="RPSFTM")
km_tse<-survfit(Surv(time=t_star,event=d_star)~trt,data=TSE$data_outcome)
OS_plot_tse<-ggsurvplot(km_tse,ylab="Overall survival",xlab=time_unit,title="TSE")

km_gform<-survfit(Surv(time=gformula_data.time/12,event=gformula_data.event)~gformula_data.trt,data=gform)
OS_plot_gform<-ggsurvplot(km_gform,ylab="Overall survival",xlab=time_unit,title="g-formula")

arrange_ggsurvplots(list(OS_plot,OS_plot_ipw,OS_plot_tse,OS_plot_gform),ncol=2,nrow=2)
#arrange_ggsurvplots(list(OS_plot_ipw,OS_plot_rps,OS_plot_tse,OS_plot_gform),ncol=2,nrow=2)

#?arrange_ggsurvplots

restab<-rbind(
  itt=itt_res,
  cens=cens_res,
  ipw=ipw_res,
  rpsftm=RPS_res,
  tse=TSE_res,
  gformula=gform_res
)
restab

#####
