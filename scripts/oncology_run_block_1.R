#Run the simulation
rm(list=ls())

sim_block<-6 #can be 1 to 6 or 0 for all
N_sim <- 10

batch_length<-8

Sim_ID<-"TEST6_"
#output path
path="results2/"

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

#Make scenario table, the code now also adds the column with scneario short names

source("scripts/oncology_make_scenario_table_1-2.R")

#Run different scenario settings:
#need to create folder "results" if not present!!
if(!file.exists("results"))  dir.create("results")

if(sim_block==1 | sim_block==10) {
  #Small n H1
  scen_select<-"small_n"
  hyp_select<-"H1"
  source("scripts/oncology_batches_1-1.R")
}

if(sim_block==2 | sim_block==0) {
  #Small n H0
  scen_select<-"small_n"
  hyp_select<-"H0"
  source("scripts/oncology_batches_1-1.R")
}

if(sim_block==3 | sim_block==0) {
  #Large n H1
  scen_select<-"large_n"
  hyp_select<-"H1"
  source("scripts/oncology_batches_1-1.R")
}

if(sim_block==4 | sim_block==0) {
  #Large n H0
  scen_select<-"large_n"
  hyp_select<-"H0"
  source("scripts/oncology_batches_1-1.R")
}

if(sim_block==5 | sim_block==0) {
  #Extra IPCW H1 (this only uses small n)
  scen_select<-"IPCW_extra"
  hyp_select<-"H1"
  source("scripts/oncology_batches_1-1.R")
}

if(sim_block==6 | sim_block==0) {
  #Extra IPCW H0 (this only uses small n)
  scen_select<-"IPCW_extra"
  hyp_select<-"H0"
  source("scripts/oncology_batches_1-1.R")
}
