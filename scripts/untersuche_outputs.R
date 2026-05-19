setwd("Z:\\EMA_causal_16May_node2\\CI.RCT.Sim\\results2")
rm(list=ls())

load("RESULTS_May_1K_scen_small_n_H1_sc_9to16_results_onco_nsim1000_ims-node32026-05-18_1922.Rdata")
library(SimDesign)
E <- SimExtract(results,"results")



load("May_1K_scen_IPCW_extra_H1_sc_1to26_results_onco_nsim1000_ims-node32026-05-18_2309.Rdata")
X <- SimExtract(results,"results")
results$rpsftm.test.rejection_0.025
results$scen_set
names(E)
sink(file="STR.txt")
str(E)
sink()
#A$´scenario´ [[1]]$rpsftm

#scenario, iteration, methode, parameter
results[[1]][[1]][["ipw"]]
X[[1]][[1]][["ipw"]]
scen<-4#12-8+1
tab<-NULL
for(i in 1:1000) {
  tab<-rbind(tab,unlist(E[[scen]][[i]][["rpsftm"]]))
}
tab<-as.data.frame(tab)

tab[is.na(tab$SElogHR),]

head(tab)
tab$p
mean(tab$p<=0.025)
is(tab)
summary(tab)

scen<-4
tabX<-NULL
for(i in 1:1000) {
  tabX<-rbind(tabX,unlist(X[[scen]][[i]][["ipw_IPCW"]]))
}
tabX<-as.data.frame(tabX)
tabX
unlist(X[[3]][[977]][["ipw_IPCW"]])
unlist(X[[4]][[977]][["ipw_IPCW"]])
unlist(X[[4]][[977]][["ipw"]])
unlist(X[[4]][[977]][["gformula"]])
unlist(X[[4]][[978]][["ipw_IPCW"]])
unlist(X[[4]][[3]][["ipw_IPCW"]])

unlist(E[[4]][[3]][["ipw"]])
unlist(E[[4]][[3]][["cens"]])
unlist(E[[4]][[3]][["cens"]])
unlist(E[[4]][[3]][["rpsftm"]])
unlist(E[[4]][[3]][["rpsftm"]])
tabX[is.na(tabX$p),]

colMeans(tabX,na.rm=TRUE)
summary(tabX)
unlist(X[[scen]][[i]][["rpsftm"]])
summary(tab)

names(results)

U<-results |>
  filter(str_detect("rpsftm", "\\.est"))
U
names(U)
