# devtools::install()
# renv::restore()
#library(CI.RCT.Sim)

rm(list=ls())
rpath<-"/home/robin/EMA_Causal_rescue_med/CI.RCT.Sim/R"
Rfiles<-list.files(path = rpath)
for(ff in Rfiles) {
  if(grepl("diabetes",ff)) source(paste(rpath,"/",ff,sep=""))

}
source(paste(rpath,"/","utility_funtions.R",sep=""))

library(future)
library(future.apply)

#library(parallel)

library(survival)
library(mvtnorm)
library(dplyr)
library(mice)
#install.packages("mice")
#library(trtswitch)
#library(SimNPH)




# Define parameter values and derived quantities -------------------------

sim_parameters <- diabetes_scenario() |>
  diabetes_scenario_set_truevalues()
sim_parameters$setup<-1
cond<-sim_parameters[1,]
cond$setup
#cond$miss[[1]][1]<- -100
data<-generate_diabetes(cond)
condition=cond
dat=data
dim(data)
summary(data)

fixe_objects=NULL
fun<-analyse_diabetes_gcomputation()
A<-fun(condition,data,NULL)
B<-fun(condition,data,NULL)


A
B

#speed

analyse_diabetes_gcomputation2 <- function() {
  function(condition, dat, fixed_objects = NULL,B=200) {
    data<-dat
    k<-condition$k
    estim_fun<-function(data) {
      mean_age<-mean(data$age)
      mean_y<-mean(data$y0)

      group<-c(0,1)
      yfinal<-rep(NA,2)
      for(u in 1:2) {
        yhat<-rep(NA,k)

        D<-data[data$trt==group[u],]

        X<-matrix(c(1,mean_y,0,mean_age),nrow=1)
        i<-1
        frml<-formula(paste("y",i,"~ y",i-1,"+R",i-1,"+age",sep=""))
        mod<-lm(frml,data=D)
        beta<-coef(mod)
        beta[is.na(beta)]<-0 #because if there is no rescue, the coefficient is NA
        yhat[1]<-X%*%beta

        for(i in 2:k) {
          X<-matrix(c(1,yhat[i-1],0,mean_y,mean_age),nrow=1)
          frml<-formula(paste("y",i,"~ y",i-1,"+R",i-1,"+y0+age",sep=""))
          mod<-lm(frml,data=D)
          beta<-coef(mod)
          beta[is.na(beta)]<-0 #because if there is no rescue, the coefficient is NA
          yhat[i]<-X%*%beta
        }
        yfinal[u]<-yhat[k]
      }
      est<-yfinal[2]-yfinal[1]
      est
    }

    est<-estim_fun(data)

    #Bootsratp
    #B<-100
    ids<-1:dim(data)[1]
    est_b<-rep(NA,B)
    for(j in 1:B) {
      ind<-sample(ids,replace=TRUE)
      data_b<-data[ind,]
      est_b[j]<-estim_fun(data_b)
    }
    SE<-sd(est_b,na.rm=TRUE)
    Q<-qt(0.975,df=B-1)
    KI<-est+c(-1,1)*Q*SE
    KI
    z<-est/SE
    p<-pnorm(z)
    p
    list(
      coef = est,
      se = SE,
      p = p,
      ci_lower = KI[1],
      ci_upper = KI[2]
    )
  }
}


fun<-analyse_diabetes_gcomputation()

fun2<-analyse_diabetes_gcomputation2()

system.time({
A<-fun(condition,data,NULL,nsamples=100)
})
system.time({
A2<-fun2(condition,data,NULL,B=100)
})
A
A2

head(data)

sim_one<-function(j) {
  data<-generate_diabetes(cond)
  A<-fun(condition,data,NULL,nsamples=10)
  A2<-fun2(condition,data,NULL,B=100)


  out<-as.data.frame(rbind(A,A2))
  rownames(out)<-c("old","new")
  out
}
sim_one(1)
start<-Sys.time()
sim_one(1)
ende<-Sys.time()
ende-start
10*1000/60/60/30




sim_one(1)
example<-sim_one(1)
example
exampleNA<-example
exampleNA[1,]<-NA
exampleNA
sim_one_try<-function(i) {
  tryCatch(sim_one(i), error = function(e) exampleNA)
}

sim_one_try(1)
plan(multisession, workers = availableCores() - 1)


R<-1000

start<-Sys.time()

res <- future_lapply(
  1:R,
  function(i) sim_one_try(i),
  future.seed = TRUE
)
end<-Sys.time()
end-start #ca. 11.2 sec for R=10


n_methods<-dim(res[[1]])[1]
n_param<-dim(res[[1]])[2]

Ar<-array(unlist(res),dim=c(n_methods,n_param,R),dimnames=list(rownames(res[[1]]),names(res[[1]]),1:R))
#T1E
Ar
dim(Ar)
apply(Ar,c(1,2),mean)
mean(Ar[1,3,]<0.025)
mean(Ar[2,3,]<0.025)

mean(Ar[1,3,]<0.1)
mean(Ar[2,3,]<0.1)


Ar[1,,]
mean(Ar[,,4])
A<-t(Ar[1,,])
A

