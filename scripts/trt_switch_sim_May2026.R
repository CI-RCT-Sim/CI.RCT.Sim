rm(list=ls())

sim_block<-1 #can be 1 to 6 or 0 for all
N_sim <- 10#00

source("scripts/oncology_make_scenario_table_1-2.R")

scen_set<-1:26 #small n, H1


Sim_ID<-"Aspera"

#
# renv::activate()
## renv::restore()
## devtools::document()
devtools::load_all()
#install.packages("fs")
#install.packages("rlang")
#install.packages("pkgload")
#devtools::install()
library(CI.RCT.Sim)
library(parallel)
library(survival)

#analysis functions
ipw_fun<-analyse_oncology_ipw2()
ipw_c_fun<-analyse_oncology_ipw2(use_censoring_IPW = TRUE, trunc_weights = 5, requ_n_cens = 5)
rpsftm_fun<-analyse_oncology_rpsftm(recensor = FALSE)
rpsftm_rc_fun<-analyse_oncology_rpsftm(recensor = TRUE)
tse_fun<-analyse_oncology_TSE(recensor = FALSE)
tse_rc_fun<-analyse_oncology_TSE(recensor = TRUE)
gformula_fun<-analyse_oncology_gformula(B = 20)
cens_fun<-analyse_oncology_cens()
rpsftm_IPCW_fun<-analyse_oncology_mixed(method="RPSFTM",recensor = TRUE,B = 50,trunc_weights = 5,use_censoring_IPW = TRUE,requ_n_cens = 5)
tse_IPCW_fun<-analyse_oncology_mixed(method="TSE",recensor = TRUE,B = 50,trunc_weights = 5,use_censoring_IPW = TRUE,requ_n_cens = 5)
gformula_IPCW_fun<-analyse_oncology_gformula(B = 20,use_censoring_IPW=TRUE, requ_n_cens=5, trunc_weights=5)
ipw_IPCW_fun<-analyse_oncology_ipw2(use_censoring_IPW = TRUE, trunc_weights = 5, requ_n_cens = 5)



make_data<-generate_oncology


sim_one_2<-function(i) {
  seed_range<-1000000
  seed<-round(runif(1,0,seed_range))
  set.seed(seed)
  data<-make_data(cond)
  res<-rbind(
    ipw=as.data.frame(ipw_fun(cond,data)),
    ipw_c=as.data.frame(ipw_c_fun(cond,data)),
    rpsftm=as.data.frame(rpsftm_fun(cond,data)),
    rpsftm_rc=as.data.frame(rpsftm_rc_fun(cond,data)),
    tse=as.data.frame(tse_fun(cond,data)),
    tse_rc=as.data.frame(tse_rc_fun(cond,data)),
    gformula=as.data.frame(gformula_fun(cond,data)),
    cens=as.data.frame(cens_fun(cond,data))#,
    #ipw_IPCW=as.data.frame(ipw_IPCW_fun(cond,data)),
    #rpsftm_IPCW=as.data.frame(rpsftm_IPCW_fun(cond,data)),
    #tse_IPCW=as.data.frame(tse_IPCW_fun(cond,data)),
    #gformula_IPCW=as.data.frame(gformula_IPCW_fun(cond,data)),
  )
  res<-as.data.frame(res)
  res$bias<-res$HR-cond$true_eff
  res$p_one_sided<-ifelse(res$HR<1,res$p/2,1-res$p/2)
  res$reject_0025<-res$p_one_sided<=0.025
  res$n_pat<-dim(data)[1]
  res$n_ctrl<-sum(data$trt==0)
  res$n_switch<-sum(data$switch)
  res$n_ev<-sum(data$ev)
  res$n_soll<-data$ev_soll[1]
  res$duration<-max(data$event_time)
  res$n_random_cens<-sum(data$random_cens)
  res$seed<-seed
  res
}
#cond=sim_parameters[1,]
#sim_one_2(1)

library(future.apply)
plan(multisession, workers = availableCores() - 4)

R<-10
R
scen_set<-1:2
scen_id<-1
for(scen_id in scen_set) {
  cond=sim_parameters[scen_id,]

start<-Sys.time()
res <- future_lapply(
  1:R,
  function(i) sim_one_2(i),
  future.seed = TRUE
)
res

end<-Sys.time()
end-start #14 min h 1000 runs with 20 cores
14000/60/24/2

n_methods<-dim(res[[1]])[1]
n_param<-dim(res[[1]])[2]

Ar<-array(unlist(res),dim=c(n_methods,n_param,R),dimnames=list(rownames(res[[1]]),names(res[[1]]),1:R))
path<-"results/"
filename=paste(path,Sim_ID,"_scen_",scen_id,"_nsim_",R,".RData",sep="")
save(Ar,file=filename)
}

#bias and coverage

dim(Ar) #Methods, Parameters, Iterations

Ar[,,1]
Ar[,,2]
resm<-apply(Ar,c(1,2),mean)
resm
Ar[,1,1]
SE_emp<-apply(log(Ar[,1,]),1,sd)

tab_sav<-cbind(resm,SE_emp,trueHR)
tab_sav

