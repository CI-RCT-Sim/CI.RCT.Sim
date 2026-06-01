#Run the simulation

#needs to be done only once
#source("scripts/oncology_get_true_values_logHR_scen_cens_L.R")
#source("scripts/oncology_get_true_values_logHR__scen_cens_L_wo_censoring.R")
rm(list=ls())

sim_block<-0 #can be 1 to 6 or 0 for all
N_sim <- 10000

batch_length<-8

Sim_ID<-"Sim_All"

#output path
path="results_sim_all"
if(!file.exists(path))  dir.create(path)


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

#test

# cl <- makeCluster(2)
# clusterEvalQ(cl, {
#   library("CI.RCT.Sim")
#   #library("survival")
# })
# alpha<-0.05
# clusterExport(cl = cl, varlist = c("alpha"))
#
# stopCluster(cl)
#

#Make scenario table, the code now also adds the column with scenario short names

trueHR_with_censoring<-read.csv("scripts/true_logHR_nsim40_n_ev10000withRandomCensTRUE.csv")
trueHR_without_censoring<-read.csv("scripts/true_logHR_nsim40_n_ev10000withRandomCensFALSE.csv")


make_scen_table<-"scripts/oncology_make_scenario_table_1-3_logHR.R"
run_batches<-"scripts/oncology_batches_1-7_logHR.R"

#save_param_tab<-TRUE
#trueHR_tab<-trueHR_with_censoring
#source(make_scen_table) #run once with saving to create the scenario list with names
save_param_tab<-FALSE

#Run different scenario settings:

if(sim_block==1 | sim_block==0) {
  #Small n H1
  scen_select<-"small_n"
  hyp_select<-"H1"
  trueHR_tab<-trueHR_with_censoring
  source(make_scen_table)
  source(run_batches)
}

if(sim_block==2 | sim_block==0) {
  #Small n H0
  scen_select<-"small_n"
  hyp_select<-"H0"
  trueHR_tab<-trueHR_with_censoring
  source(make_scen_table)
  source(run_batches)
}

if(sim_block==3 | sim_block==0) {
  #Large n H1
  scen_select<-"large_n"
  hyp_select<-"H1"
  trueHR_tab<-trueHR_with_censoring
  source(make_scen_table)
  source(run_batches)
}

if(sim_block==4 | sim_block==0) {
  #Large n H0
  scen_select<-"large_n"
  hyp_select<-"H0"
  trueHR_tab<-trueHR_with_censoring
  source(make_scen_table)
  source(run_batches)
}

if(sim_block==5 | sim_block==0) {
  #Extra IPCW H1 (this only uses small n)
  scen_select<-"IPCW_extra"
  hyp_select<-"H1"
  trueHR_tab<-trueHR_without_censoring
  source(make_scen_table)
  source(run_batches)
}

if(sim_block==6 | sim_block==0) {
  #Extra IPCW H0 (this only uses small n)
  scen_select<-"IPCW_extra"
  hyp_select<-"H0"
  trueHR_tab<-trueHR_without_censoring
  source(make_scen_table)
  source(run_batches)
}
