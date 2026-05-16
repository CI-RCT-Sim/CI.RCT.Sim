
#to be run from oncology_run_all.R

#devtools::load_all()
#rm(list=ls())
#library(CI.RCT.Sim)
#library(parallel)

pre_sim_parameters <- oncology_scenario()
A<-pre_sim_parameters |> oncology_scenario_set_truevalues()
head(A)

A$beta_cens

#beta_lab<-c("Int","X","W","Wgrw","L","trt","switched")
#beta_lab2<-c("Int","X","W","Wgrw","L","trt","switched","logHR_assumed")
#beta_lab3<-c("Int","X","W","Wgrw","L","trt","switched","switching_in_control_only")


H0<-rep(FALSE,dim(A)[1])
for(i in 1:length(H0)) {
 H0[i]<-A$beta_death[[i]][6]==0
}
H0
small_n<-A$ev_soll==66

all_param_tab<-NULL

LIST<-list(
  scen_trt_sw_H1_high=A[!H0 & small_n,],
  scen_trt_sw_H1_low=A[!H0 & !small_n,],
  scen_trt_sw_H0_high=A[H0 & small_n,],
  scen_trt_sw_H0_low=A[H0 & !small_n,]
)
block_nam<-names(LIST)[1]
for(block_nam in names(LIST)) {
  tab<-LIST[[block_nam]]

  #mu W:
  temp<-tab$mu_W
  j<-1
  pat<-NULL
  for(j in 1:length(temp)) {
    pat<-c(pat,paste(temp[[j]]$trt,collapse=", "))
    pat<-c(pat,paste(temp[[j]]$ctr,collapse=", "))
  }
  pat_num<-as.numeric(factor(pat))
  muW_levels<-levels(factor(pat))
  muW_pattern<-matrix(pat_num,ncol=2,byrow=TRUE)
  colnames(muW_pattern)<-c("W_mean_pattern_Trt","W_mean_pattern_Ctr")
  muW_pattern

  #mu L
  temp<-tab$mu_L
  j<-1
  pat<-NULL
  for(j in 1:length(temp)) {
    pat<-c(pat,paste(temp[[j]]$trt,collapse=", "))
    pat<-c(pat,paste(temp[[j]]$ctr,collapse=", "))
  }
  pat_num<-as.numeric(factor(pat))
  muL_levels<-levels(factor(pat))
  muL_pattern<-matrix(pat_num,ncol=2,byrow=TRUE)
  colnames(muL_pattern)<-c("L_mean_pattern_Trt","L_mean_pattern_Ctr")
  muL_pattern

  #W and L covariance matrix
  temp<-tab$Sigma_W_L
  j<-1
  patVmat<-NULL
  #Pat_Sigma<-NULL
  for(j in 1:length(temp)) {
    #Pat_Sigma<-rbind(Pat_Sigma,temp[[j]])
    patVmat<-c(patVmat,paste(temp[[j]],collapse=", "))
  }
  patVmat_num<-as.numeric(factor(patVmat))
  Sigma_levels<-levels(factor(patVmat))
  Vmat_pattern<-matrix(patVmat_num,ncol=1,byrow=TRUE)
  colnames(Vmat_pattern)<-c("Covariance matrix type")
  Vmat_pattern


  #Betas:

  #beta_lab<-c("Int","X","W","Wgrw","L","trt","switched")
  #beta_lab2<-c("Int","X","W","Wgrw","L","trt","switched","logHR_assumed")
  #beta_lab3<-c("Int","X","W","Wgrw","L","trt","switched","switching_in_control_only")

  BETA<-NULL
  beta_names<-names(tab)[grepl("beta",names(tab))]
  b<-beta_names[1]
  for(b in beta_names) {
    temp<-tab[[b]]
    betas<-matrix(unlist(temp),ncol=length(temp[[1]]),byrow=TRUE)
    colnames(betas)<-paste(b,names(temp[[1]]),sep=".")
    #if(b=="beta_prog" | b=="beta_switch") colnames(betas)<-beta_lab
    #if(b=="beta_death") colnames(betas)<-beta_lab2
    #if(b=="beta_cens") colnames(betas)<-beta_lab3
    #colnames(betas)<-paste(b,colnames(betas),sep=".")
    BETA<-cbind(BETA,betas)
  }
  BETA<-exp(BETA) #to have Hazard ratios in table
  #constants:
  names(tab)
  constant<-tab[,c("k","recr_interval","max_duration","aimed_for_n_per_required_event","alpha","power","p_trt","w")]

  param_tab<-data.frame(block_nam=block_nam,muW_pattern,muL_pattern,Vmat_pattern,BETA,constant)

  all_param_tab<-rbind(all_param_tab,param_tab)
}

all_param_tab<-cbind(Scen_ID=1:dim(all_param_tab)[1],all_param_tab)



names(all_param_tab)
all_param_tab$beta_cens.switching_in_control_only<-round(log(all_param_tab$beta_cens.switching_in_control_only))

#names:
names_temp<-c(
  "Core",
  "W - decrease under ctr",
  "W - decrease in both groups",
  "L - decrease under ctr, effect on death, progr., switching (unobs. conf.)",
  "L - decrease in both groups, effect on death, progr., switching (unobs. conf.)",
  "time dependent correlation - Toeplitz",
  "progression - fast",
  "progresssion - no effect of W",
  "progression - reduced effect of trt",
  "progression - no effect of trt",
  "progression - unobserved confounding with death",
  "switch - high probability",
  "switch - no effect of W",
  "switch - W>0 near separation",
  "switch - unobserved confounding with death",
  "death - high rate",
  "death - no effect of W",
  "death - unobserved confounding with progr. and switch",
  "Trt effect post switch - reduced",
  "Trt effect post switch - none",
  "random censoring - none",
  "random censoring - effect of X",
  "random censoring - effect of X and W",
  "random censoring - effect of X, W and L",
  "random censoring - high rate, effect of X, W",
  "random censoring - high rate, control only, effect of X, W")

scenario_names<-c(
  names_temp,
  names_temp[names_temp!="Trt effect post switch - reduced"],
  names_temp[!grepl("Trt effect post switch",names_temp)],
  names_temp[!grepl("Trt effect post switch",names_temp)]
)


all_param_tab$scenario_name<-scenario_names
dim2<-dim(all_param_tab)[2]
all_param_tab<-all_param_tab[,c(1,2,dim2,3:(dim2-1))]
head(all_param_tab)
#scen_tab<-all_param_tab

#true values were calculated using oncology_get_true_values.R, note this uses scneario tables 1-1,
#which includes the H0 scnario with unequal trajectories for W
#but this scenario is removed at the end of this code, because it is not really H0
trueHR<-read.csv("scripts/trueHR_nsim20_n_ev10000.csv")
trueHR<-trueHR[,c("scen","trueHR")]
names(trueHR)<-c("Scen_ID","true_eff")
head(trueHR)
dim(trueHR)
dim(all_param_tab)
scen_tab<-merge(all_param_tab,trueHR,by="Scen_ID",all.x=TRUE,all.y=FALSE)
head(scen_tab)
dim(scen_tab)
scen_tab$true_eff[is.na(scen_tab$true_eff)]<-1

#remove H0 scenario with unequal trajectories for W
head(scen_tab)
rem<-grepl("H0",scen_tab$block_nam) & grepl("W - decrease under ctr",scen_tab$scenario_name)
rem
scen_tab<-scen_tab[!rem,]
#dim(scen_tab)
