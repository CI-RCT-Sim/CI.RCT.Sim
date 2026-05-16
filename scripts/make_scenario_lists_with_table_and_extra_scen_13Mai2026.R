#rm(list=ls())


params_scenarios_grid <- function(...){
  params_args <- list(...)
  params_ref <- purrr::map(params_args, \(x){
    x[1]
  }) |>
    tibble::as_tibble()

  params_other <- purrr::imap(params_args, \(x, i){
    if(inherits(x, "list")){
      purrr::map(x[-1], \(y){
        tmp <- params_ref
        tmp[,i] <- list(list(y))
        tmp
    }) |>
      purrr::list_rbind()
    } else {
      purrr::map(x[-1], \(y){
        tmp <- params_ref
        tmp[,i] <- y
        tmp
      }) |>
        purrr::list_rbind()
    }
  }) |>
    purrr::list_rbind()

  rbind(params_ref, params_other)
}




AR_1_fun<-function(rho,k) {
  i<-rep(1:k,k)
  j<-rep(1:k,each=k)
  d<-abs(i-j)
  corr<-rho^d
  matrix(corr,k,k)
}
#AR_1_fun(0.9,10)

exch_fun<-function(rho,k) {
  M<-matrix(rho,k,k)
  diag(M)<-1
  M
}

toeplitz_fun<-function(decr=0.1,k) {
  korr_vek<-seq(from=1, by=-decr, length.out=k)
  korr_vek[korr_vek<0]<-0
  i<-rep(1:k,k)
  j<-rep(1:k,each=k)
  d<-abs(i-j)+1
  corr<- korr_vek[d]
  matrix(corr,k,k)
}
#toeplitz_fun(0.1,10)


k<-10

label_fun<-function(x,labels) {
	for(i in 1:length(x)) names(x[[i]])<-labels
	x
}



mu_W<-list(
	list(trt=rep(0,k),ctr=rep(0,k)),
	list(trt=rep(0,k),ctr=c(1,0.5,0,rep(-1,k-3))),
	list(trt=c(1,0.5,0,rep(-1,k-3)),ctr=c(1,0.5,0,rep(-1,k-3)))
)
mu_L<-mu_W
Sigma_W_L<-list(exch_fun(0.5,k),toeplitz_fun(0.1,k))

beta_lab<-c("Int","X","W","Wgrw","L","trt","switched")
beta_lab2<-c("Int","X","W","Wgrw","L","trt","switched","logHR_assumed")
beta_lab3<-c("Int","X","W","Wgrw","L","trt","switched","swtiching_in_control_only")


beta_prog<-list(
	#    Int,             X,          W,    W>0, L,       trt ,      switched
	c( log( log(2)/0.5), log(0.5), log(0.5),  0, 0,        log(0.5),  0 ),
	c( log( log(2)/1)  , log(0.5), log(0.5),  0, 0,        log(0.5),  0 ),
	c( log( log(2)/0.5), log(0.5), 0       ,  0, 0,        log(0.5),  0 ),
	c( log( log(2)/0.5), log(0.5), log(0.5),  0, log(0.5), log(0.5),  0 ),
	c( log( log(2)/0.5), log(0.5), log(0.5),  0, log(0.5), log(0.75), 0 ),
	c( log( log(2)/0.5), log(0.5), log(0.5),  0, log(0.5), 0        , 0 )
) |> label_fun(beta_lab)

beta_switch<-list(
	#   Int,          X,       W,      W>0,                      L   ,trt,  ,switched
	c(log(0.5/0.5),log(1.5),log(1.5),0                         ,0       ,0,0),
	c(log(0.9/0.1),log(1.5),log(1.5),0                         ,0       ,0,0),
	c(log(0.5/0.5),log(1.5),0       ,0                         ,0       ,0,0),
	c(log(0.5/0.5),log(1.5),log(1.5),log(0.9/0.1)*pi/2-log(1.5),0       ,0,0),
	c(log(0.5/0.5),log(1.5),log(1.5),0                         ,log(1.5),0,0)
) |> label_fun(beta_lab)

HR_h<-0.5 #HR assumed high effect

beta_death_high_effect<-list(
	#     Int,          X,          W,    W>0,  L,          trt,        switched,   HR_assumed
	c( log( log(2)/4), log(0.5), log(0.5),  0, 0       ,  log(0.5) ,    log(0.5) , log(HR_h) ),
	c( log( log(2)/2), log(0.5), log(0.5),  0, 0       ,  log(0.5) ,    log(0.5) , log(HR_h) ),
	c( log( log(2)/4), log(0.5), 0       ,  0, 0       ,  log(0.5) ,    log(0.5) , log(HR_h) ),
	c( log( log(2)/4), log(0.5), log(0.5),  0, log(0.5),  log(0.5) ,    log(0.5) , log(HR_h) ),
	c( log( log(2)/4), log(0.5), log(0.5),  0, 0       ,  log(0.5) ,    0        , log(HR_h) ),
	c( log( log(2)/4), log(0.5), log(0.5),  0, 0       ,  log(0.5) ,    log(0.75), log(HR_h) )
) |> label_fun(beta_lab2)


HR_l<-0.75 #HR assumed low effect
beta_death_low_effect<-list(
	#     Int,          X,          W,    W>0,  L,          trt,        switched   ,   HR_assumed
	c( log( log(2)/4), log(0.5), log(0.5),  0, 0       ,  log(0.75) ,    log(0.75) , log(HR_l)  ),
	c( log( log(2)/2), log(0.5), log(0.5),  0, 0       ,  log(0.75) ,    log(0.75) , log(HR_l)  ),
	c( log( log(2)/4), log(0.5), 0       ,  0, 0       ,  log(0.75) ,    log(0.75) , log(HR_l)  ),
	c( log( log(2)/4), log(0.5), log(0.5),  0, log(0.5),  log(0.75) ,    log(0.75) , log(HR_l)  ),
	c( log( log(2)/4), log(0.5), log(0.5),  0, 0       ,  log(0.75) ,    0         , log(HR_l)  )
) |> label_fun(beta_lab2)



beta_death_H0_high_effect<-list(
	#    Int,              X,       W,    W>0,   L,       trt, switched   , HR_assumed
	c( log( log(2)/4), log(0.5), log(0.5),  0, 0       ,  0 ,    0 , log(HR_h)  ),
	c( log( log(2)/2), log(0.5), log(0.5),  0, 0       ,  0 ,    0 , log(HR_h)  ),
	c( log( log(2)/4), log(0.5), 0       ,  0, 0       ,  0 ,    0 , log(HR_h)  ),
	c( log( log(2)/4), log(0.5), log(0.5),  0, log(0.5),  0 ,    0 , log(HR_h)  )
) |> label_fun(beta_lab2)



beta_death_H0_low_effect<-list(
	#    Int,              X,       W,     W>0,  L,       trt, switched   , HR_assumed
	c( log( log(2)/4), log(0.5), log(0.5),  0, 0       ,  0 ,    0 , log(HR_l)  ),
	c( log( log(2)/2), log(0.5), log(0.5),  0, 0       ,  0 ,    0 , log(HR_l)  ),
	c( log( log(2)/4), log(0.5), 0       ,  0, 0       ,  0 ,    0 , log(HR_l)  ),
	c( log( log(2)/4), log(0.5), log(0.5),  0, log(0.5),  0 ,    0 , log(HR_l)  )
)  |> label_fun(beta_lab2)


beta_cens<-list(
	#    Int,                X,       W,       W>0, L,        trt, switched, switching in control only
	c(log(-log(1-0.025))  ,0       , 0       ,  0,  0       ,    0,  0,		0),
	c(-Inf                ,0       , 0       ,  0,  0       ,    0,  0,		0),
	c(log(-log(1-0.025))  ,log(0.5), 0       ,  0,  0       ,    0,  0,		0),
	c(log(-log(1-0.025))  ,log(0.5), log(0.5),  0,  0       ,    0,  0,		0),
	c(log(-log(1-0.025))  ,log(0.5), log(0.5),  0,  log(0.5),    0,  0,		0),
	#extra
	c(log(-log(1-0.1))  ,log(0.5), log(0.5),  0,  0       ,    0,  0,			0),
	c(log(-log(1-0.1))  ,log(0.5), log(0.5),  0,  0       ,    0,  0,			1)


) |> label_fun(beta_lab3)



constant<-data.frame(
	k=k,
	recr_interval=2,
	max_duration=7,
	alpha=0.05,
	power=0.8,
	p_trt=0.5,
	w=0,
	aimed_for_n_per_required_event=3
)



scen_trt_sw_H1_high<-cbind(params_scenarios_grid(mu_W=mu_W,mu_L=mu_L,Sigma_W_L=Sigma_W_L,beta_prog=beta_prog,beta_switch=beta_switch,beta_death=beta_death_high_effect,beta_cens=beta_cens),constant)
scen_trt_sw_H1_low<-cbind(params_scenarios_grid(mu_W=mu_W,mu_L=mu_L,Sigma_W_L=Sigma_W_L,beta_prog=beta_prog,beta_switch=beta_switch,beta_death=beta_death_low_effect,beta_cens=beta_cens),constant)
scen_trt_sw_H0_high<-cbind(params_scenarios_grid(mu_W=mu_W,mu_L=mu_L,Sigma_W_L=Sigma_W_L,beta_prog=beta_prog,beta_switch=beta_switch,beta_death=beta_death_H0_high_effect,beta_cens=beta_cens),constant)
scen_trt_sw_H0_low<-cbind(params_scenarios_grid(mu_W=mu_W,mu_L=mu_L,Sigma_W_L=Sigma_W_L,beta_prog=beta_prog,beta_switch=beta_switch,beta_death=beta_death_H0_low_effect,beta_cens=beta_cens),constant)

all_scen<-rbind(
  scen_trt_sw_H1_high,
  scen_trt_sw_H1_low,
  scen_trt_sw_H0_high,
  scen_trt_sw_H0_low
)
#ev_soll<-ceiling(((qnorm(1-param$alpha/2) + qnorm(param$power)) / logHR_assumed )^2 /param$p_trt/(1-param$p_trt))



dim(scen_trt_sw_H1_high)
dim(scen_trt_sw_H1_low)
dim(scen_trt_sw_H0_high)
dim(scen_trt_sw_H0_low)
length(scen_trt_sw_H1_high)

#str(scen_trt_sw_H1_high[1,])
#param<-scen_trt_sw_H1_high[1,]
#str(param)
#is(param$beta_prog[[1]])
#param$beta_prog[[1]]["X"]
#param$beta_cens[[1]]["Int"]<-log(-log(1-0.025))


#names(scen_trt_sw_H0_high[9,])

#scen_trt_sw_H0_high[9,]$beta_prog
#scen_trt_sw_H1_high[1,]$beta_prog

#scen_trt_sw_H0_high[9,]$beta_switch
#scen_trt_sw_H1_high[1,]$beta_switch


str(scen_trt_sw_H1_high)
setwd("Z:\\Projekte\\EMA_Tender_causal\\temp")
scen_trt_sw_H1_high[[1]]

###

all_param_tab<-NULL

LIST<-list(
	scen_trt_sw_H1_high=scen_trt_sw_H1_high,
	scen_trt_sw_H1_low=scen_trt_sw_H1_low,
	scen_trt_sw_H0_high=scen_trt_sw_H0_high,
	scen_trt_sw_H0_low=scen_trt_sw_H0_low
)
block_nam<-names(LIST)[1]
for(block_nam in names(LIST)) {
tab<-LIST[[block_nam]]

#mu W:
temp<-tab[[1]]
j<-1
pat<-NULL
for(j in 1:length(temp)) {
	pat<-c(pat,paste(temp[[j]]$trt,collapse=", "))
	pat<-c(pat,paste(temp[[j]]$ctr,collapse=", "))
}
pat_num<-as.numeric(factor(pat))
levels(factor(pat))
muW_pattern<-matrix(pat_num,ncol=2,byrow=TRUE)
colnames(muW_pattern)<-c("W_mean_pattern_Trt","W_mean_pattern_Ctr")
muW_pattern

#mu L
temp<-tab[[2]]
j<-1
pat<-NULL
for(j in 1:length(temp)) {
	pat<-c(pat,paste(temp[[j]]$trt,collapse=", "))
	pat<-c(pat,paste(temp[[j]]$ctr,collapse=", "))
}
pat_num<-as.numeric(factor(pat))
levels(factor(pat))
muL_pattern<-matrix(pat_num,ncol=2,byrow=TRUE)
colnames(muL_pattern)<-c("L_mean_pattern_Trt","L_mean_pattern_Ctr")
muL_pattern

#W and L covariance matrix
temp<-tab[[3]]
j<-1
patVmat<-NULL
for(j in 1:length(temp)) {
	patVmat<-c(patVmat,paste(temp[[j]],collapse=", "))
}
patVmat_num<-as.numeric(factor(patVmat))
levels(factor(patVmat))
Vmat_pattern<-matrix(patVmat_num,ncol=1,byrow=TRUE)
colnames(Vmat_pattern)<-c("Covariance matrix type")
Vmat_pattern


#Betas:
BETA<-NULL
beta_names<-names(tab)[grepl("beta",names(tab))]
b<-beta_names[1]
for(b in beta_names) {
	temp<-tab[[b]]
	betas<-matrix(unlist(temp),ncol=length(temp[[1]]),byrow=TRUE)
	colnames(betas)<-paste(b,names(temp[[1]]),sep=".")
	BETA<-cbind(BETA,betas)
}
BETA<-exp(BETA) #to have Hazard ratios in table
#constants:
param_tab<-data.frame(block_nam=block_nam,muW_pattern,muL_pattern,Vmat_pattern,BETA,constant)

all_param_tab<-rbind(all_param_tab,param_tab)
}

all_param_tab<-cbind(Scen_ID=1:dim(all_param_tab)[1],all_param_tab)

#library(writexl)
setwd("Z:\\Projekte\\EMA_Tender_causal\\Simulation")
#write_xlsx(list(full_data=tab2,nur_jene_mit_VARK=tab_restr2),"Gruppenvergleich_mean_SD.xlsx")
#write_xlsx(all_param_tab,"all_param_tab_temp.xlsx")

for(i in 2:25) {
print(which(all_param_tab[i,]!=all_param_tab[1,]))
}


for(i in 27:49) {
print(which(all_param_tab[i,]!=all_param_tab[26,]))
}


for(i in 51:72) {
print(which(all_param_tab[i,]!=all_param_tab[50,]))
}



which(all_param_tab[4,]!=all_param_tab[1,])
which(all_param_tab[5,]!=all_param_tab[1,])
which(all_param_tab[7,]!=all_param_tab[1,])




library(openxlsx)

# write dataset
wb <- createWorkbook()
addWorksheet(wb, sheetName="Scenarios")
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
getwd()
saveWorkbook(wb, "yellow_7May.xlsx", overwrite=TRUE)
