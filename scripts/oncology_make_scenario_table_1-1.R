
devtools::load_all()
rm(list=ls())
library(CI.RCT.Sim)
library(parallel)

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



library(openxlsx)

# write dataset
wb <- createWorkbook()
addWorksheet(wb, sheetName="Scenarios")
names(all_param_tab)
all_param_tab$beta_cens.switching_in_control_only<-round(log(all_param_tab$beta_cens.switching_in_control_only))

writeData(wb, sheet="Scenarios", x=all_param_tab)

# define style
yellow_style <- createStyle(fgFill="#FFFF00")

# difference to core scenario:
checko<-matrix(FALSE,nrow=dim(all_param_tab)[1],ncol=dim(all_param_tab)[2])
b<-names(LIST)[1]

for(b in names(LIST)) {
  set<-all_param_tab$block_nam==b
  ind<-(1:dim(all_param_tab)[1])[set]
  for(i in ind[-1]) {
    different<-which(all_param_tab[i,]!=all_param_tab[ind[1],])
    checko[i,different]<-TRUE
  }
}



for(x in 1:dim(all_param_tab)[1]) {
  for(y in 2:dim(all_param_tab)[2]) { #start at 2, because 1 is the ID and these are all different but should not be highlighted
    if(checko[x,y]) addStyle(wb, sheet="Scenarios", style=yellow_style, rows=x+1, cols=y, gridExpand=TRUE) # +1 for header line
  }
}
# write result


#
addWorksheet(wb, sheetName="Mean_W_patterns")
Pat_muW<-data.frame(Pattern=1:length(muW_levels),Value=muW_levels)
writeData(wb, sheet="Mean_W_patterns", x=Pat_muW)

addWorksheet(wb, sheetName="Mean_L_patterns")
Pat_muL<-data.frame(Pattern=1:length(muL_levels),Value=muL_levels)
writeData(wb, sheet="Mean_L_patterns", x=Pat_muL)

addWorksheet(wb, sheetName="Covariance_patterns")
matrix(unlist(Sigma_levels),ncol=sqrt(length(Sigma_levels[[1]])))

Pat_Sigma<-NULL
for(i in 1:length(Sigma_levels)) {
  xx<-Sigma_levels[i]
  yy<-as.numeric(strsplit(xx,", ")[[1]])
  mat<-matrix(yy,ncol=sqrt(length(yy)))
  colnames(mat)<-paste("V",1:dim(mat)[2],sep="")
  Pat_Sigma<-rbind(Pat_Sigma,cbind(Pattern=i,mat))
}

#Pat_Sigma
writeData(wb, sheet="Covariance_patterns", x=Pat_Sigma)

saveWorkbook(wb, "yellow_13May_1-3.xlsx", overwrite=TRUE)













#
