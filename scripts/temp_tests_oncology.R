#Test
rm(list=ls())

devtools::load_all()
library(CI.RCT.Sim)
library(parallel)

sim_parameters <- oncology_scenario()
is(sim_parameters)
A<-sim_parameters[27,]
#A<-as.data.frame(A)
data<-generate_oncology(A)
head(data)


analyse_oncology_ipw()(A,data)

system.time(analyse_oncology_gformula(B=20)(A,data))
system.time(analyse_oncology_gformula(B=200)(A,data))

analyse_oncology_gformula(B=20,n_ev_cutoff_no_bootstrap=100)(A,data)



B<-all_scen[90,]

B

data<-generate_oncology(B)

analyse_oncology_ipw()(para,data)

str(A)
str(B)

sort(names(A))
sort(names(B))
A$aimed_for_n_per_required_event
