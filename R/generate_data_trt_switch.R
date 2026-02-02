#rm(list=ls())
#treatment_switching_1-11.R

library(mvtnorm)

expit<-function(x) exp(x)/(1+exp(x))


AR_1_fun<-function(rho,k) {
  i<-rep(1:k,k)
  j<-rep(1:k,each=k)
  d<-abs(i-j)
  corr<-rho^d
  matrix(corr,k,k)
}
exch_fun<-function(rho,k) {
  M<-matrix(rho,k,k)
  diag(M)<-1
  M
}


rtime<-function(n,TRT,SWITCHED,X,W,L,b0,b_trt,b_sw,b_x,b_W,b_L) {

  #event times
  zz<-rep(NA,n)
  u<- -log(runif(n))  #Inverse CDF transformed to cumulative hazard
  for(ii in 1:n) {
    eta<-exp(b0+TRT[ii,]*b_trt+SWITCHED[ii,]*b_sw+X[ii,]*b_x+W[ii,]*b_W+L[ii,]*b_L)
    #the term SWITCHED[ii,]*b_sw only provides an approximation, need to refine
    cumhaz<-cumsum(eta) #cumulative hazard after 1, 2, 3, ..., k years. Values before, between or after these timepoints are interplated in the next steps. Note that the
    pos<-sum(cumhaz<=u[ii])
    if(pos==0) eta_pos<-0 else eta_pos<-sum(eta[1:pos])
    if(pos==k) eta<-c(eta,eta[k]) #so to define eta[pos+1]
    zz[ii] <- (u[ii] - eta_pos)/eta[pos+1] + pos
  }
  zz
}

#Parameters
k<-10
rho_obs<-0.9 #AR(1)
rho_unobs<-0.9 #AR(1)
sd_obs<-1
sd_unobs<-1


param_list_trt_sw<-list(scenario_1=
                          list(

                            #maximal number of visits
                            k=k,
                            #rho=0.9,
                            #Study time frames:
                            recr_interval=2,
                            max_duration=6,
                            #baseline hazard
                            b0=rep( log( log(2)/4  ),k), #median 4
                            g0=rep( log( log(2)/0.5  ),k), #median 1
                            #Treatment effect (log HR),
                            b_trt= log(0.5),
                            g_trt= log(0.5),

                            b_x=0.5, #b...coefficients (log HR) for hazard of death
                            g_x=0.5, #g...coefficients (log HR) for hazard of progression
                            #baseline covariate x
                            mu_x=0,
                            sd_x=1,

                            #observed time dependent covariates W
                            mu_obs=rep(0,k),
                            V_obs=AR_1_fun(rho_obs,k)*sd_obs^2,
                            b_W= rep(0.5,k),
                            g_W=rep(0.5,k),

                            #unobserved  time dependent covariates L
                            mu_unobs=rep(0,k),
                            V_unobs=AR_1_fun(rho_unobs,k)*sd_unobs^2,
                            b_L=rep(0.5,k),
                            g_L=rep(0.5,k),

                            #Effect of active treatment after switching
                            g_sw=rep(0,k), #does not matter, because for the progression mechanism, switching is not used, as switching can only occur after progression
                            b_sw=log(0.5), #as b_trt or smaller

                            #switching
                            p_switch_0=0.5,
                            h_0=0, #log(p_switch_0/(1-p_switch_0))
                            h_X=log(1.5),
                            h_W=log(1.5),

                            alpha=0.05,
                            power=0.8,
                            p_trt=0.5,
                            HR_assumed=0.5  #to determine sample size
                          )
)


#Scenario 1
make_data_trt_switch<-function(param,allow_switch=TRUE) {
  k=param$k
  recr_interval=param$recr_interval
  max_duration=param$max_duration
  b0=param$b0
  g0=param$g0
  b_trt=param$b_trt
  g_trt=param$g_trt
  b_x=param$b_x
  g_x=param$g_x
  mu_x=param$mu_x
  sd_x=param$sd_x
  mu_obs=param$mu_obs
  V_obs=param$V_obs
  b_W=param$b_W
  g_W=param$g_W
  mu_unobs=param$mu_unobs
  V_unobs=param$V_unobs
  b_L=param$b_L
  g_L=param$g_L
  g_sw=param$g_sw
  b_sw=param$b_sw
  p_switch_0=param$p_switch_0
  h_0=param$h_0
  h_X=param$h_X
  h_W=param$h_W
  alpha=param$alpha
  power=param$power
  p_trt=param$p_trt
  HR_assumed=param$HR_assumed

  #
  ev_soll<-((qnorm(1-alpha/2) + qnorm(power)) / log(HR_assumed) )^2 /p_trt/(1-p_trt)
  n_aim<-ceiling(ev_soll)*2
  n<-rpois(1,n_aim)

  X<-matrix(rep(rnorm(n,mu_x,sd_x),k),ncol=k)
  W<-rmvnorm(n,mean=mu_obs,sigma=V_obs)
  L<-rmvnorm(n,mean=mu_unobs,sigma=V_unobs)

  trt<-rbinom(n,size=1,prob=p_trt)
  TRT<-matrix(rep(trt,k),ncol=k)
  SWITCHED<-matrix(0,nrow=n,ncol=k)

  #SWITCHED oder TRT muss spater fractions erlauben, fuer das Intervall, in dem geswitched wird. Move to analysis function
  prog_time<-rtime(n,TRT,SWITCHED,X,W,L,g0,g_trt,g_sw,g_x,g_W,g_L)

  #covariate values at progression (secondary baseline)
  index_sec_BL<-floor(prog_time)+1
  index_sec_BL[index_sec_BL>k ]<-k

  X_2BL<-W_2BL<-L_2BL<-TRT_2BL<-rep(NA,n)
  for(j in 1:n) {
    X_2BL[j]<-X[j,index_sec_BL[j]]
    W_2BL[j]<-W[j,index_sec_BL[j]]
    L_2BL[j]<-L[j,index_sec_BL[j]]
  }


  if(allow_switch) {
    switch_prob<-ifelse(trt==0,expit(h_0 + X_2BL*h_X + W_2BL*h_W),0)
  } else {
    switch_prob<-rep(0,n)
  }
  switch<-rbinom(n,size=1,prob=switch_prob)==1
  for(j in (1:n)[switch & index_sec_BL<k]) {
    #TRT[j,(index_sec_BL[j]+1):k]<-1
    #TRT[j,index_sec_BL[j]]<- prog_time[j]%%1  #fraction of time under trt=1 in this interval
    SWITCHED[j,(index_sec_BL[j]+1):k]<-1
    SWITCHED[j,index_sec_BL[j]]<- prog_time[j]%%1
  }

  for(j in 1:n) {
    TRT_2BL[j]<-TRT[j,index_sec_BL[j]] + SWITCHED[j,index_sec_BL[j]]
  }

  ind_PD<-0:(k-1)
  PD<-outer(prog_time,ind_PD,"<=")

  ind<-0:(k-1) #used for variable names to indicate visits. Visit 0 is baseline, then visits are performed every year. So ind is the time of the visit in years.

  event_time_uncensored<-rtime(n,TRT,SWITCHED,X,W,L,b0,b_trt,b_sw,b_x,b_W,b_L)

  start<-runif(n,0,recr_interval)
  cal<-start+event_time_uncensored #calendar times
  colnames(TRT)<-paste("TRT",ind,sep="_")
  colnames(X)<-paste("X",ind,sep="_")
  colnames(W)<-paste("W",ind,sep="_")
  colnames(L)<-paste("L",ind,sep="_")
  #colnames(switch)<-paste("switch",ind,sep="_")
  PD<-matrix(as.numeric(PD),nrow=dim(PD)[1],ncol=dim(PD)[2])
  colnames(PD)<-paste("PD",ind,sep="_")
  colnames(SWITCHED)<-paste("SWITCHED",ind,sep="_")

  temp<-data.frame(
    id=1:n,
    trt,
    #TRT,
    X,W,L,
    #PD,switch,SWITCHED,index_sec_BL,TRT_2BL,X_2BL,W_2BL,L_2BL,
    cal,
    event_time_uncensored,
    start,
    event_time=NA,
    ev=1,
    prog_time,
    prog_ev=NA,
    switch
  )

  #administrative censoring
  temp<-temp[order(temp$cal),]
  #temp$ev<-1
  temp$cum_ev<-cumsum(temp$ev)
  cal_end<-temp$cal[ev_soll]

  temp$ev<-as.numeric(temp$cal<=cal_end)
  temp$cal2<-ifelse(temp$ev==1,temp$cal,cal_end) #censored calendar time

  temp$event_time<-temp$cal2-temp$start

  #temp$sec_BL_reached<-temp$y>(temp$index_sec_BL-1)

  #temp$prog_time<=temp$event_time
  temp$prog_ev<-as.numeric(temp$prog_time<temp$event_time)
  temp$prog_time[!temp$prog_ev]<-temp$event_time[!temp$prog_ev]

  #temp$PFS_event<-temp$prog_ev | temp$ev
  #temp$PFS_time<-pmin(temp$prog_time,temp$event_time)

  temp$calendar_start_time<-temp$start
  temp$calendar_end_of_study<-cal_end
  #temp$max_FU<-cal_end-temp$start #the maximal follow up to be used for RPSFTM recensoring, but we should calculate this in the analysis functions


  temp$cal<-NULL
  temp$event_time_uncensored<-NULL
  temp$start<-NULL
  temp$cum_ev<-NULL
  temp$cal2<-NULL

  temp$switch<-as.numeric(temp$switch & temp$prog_ev==1)

  #remove covariate values at unobserved timepoints
  ind<-0:(k-1)
  i<-1
  for(i in 1:dim(temp)[1]) {
    unobs<-ind>temp$event_time[i]
    rm_X<-paste(rep(c("X","W","L"),each=sum(unobs)),ind[unobs],sep="_")
    temp[i,rm_X]<-NA
  }

  temp
}

#Example use
if(FALSE) {
  data<-make_data_trt_switch(param_list_trt_sw$scenario_1)
  head(data)
  library(survival)
  library(survminer)

  km_OS<-survfit(Surv(time=event_time,event=ev)~trt,data=data)
  OS_plot<-ggsurvplot(km_OS,data=data,ylab="Overall survival",xlab="Years",risk.table=TRUE)
  OS_plot

  coxph(Surv(time=event_time,event=ev)~trt+X_0+W_0,data=data)
}



