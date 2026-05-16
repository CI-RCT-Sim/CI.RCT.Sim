#Run the simulation

# devtools::install()
# renv::activate()
## renv::restore()
## devtools::document()
devtools::load_all()
rm(list=ls())
library(CI.RCT.Sim)
library(parallel)
library(survival)

#Make scenario table, the code now also adds the column with scneario short names
#save_param_tab<-FALSE #need to create folder "data" if not present, as it is saved there
source("scripts/oncology_make_scenario_table_1-2.R")

#Run different scenario settings:
#need to create folder "results" if not present!!
if(!file.exists("results"))  dir.create("results")


# Iterations and scenarios
N_sim <- 10000 #maybe try it first with a small number like 10


#### run with these three settings separately:

batch_length<-8

#Small n H1
scen_select<-"small_n"
hyp_select<-"H1"
source("scripts/oncology_batches_1-1.R")

#Small n H0
scen_select<-"small_n"
hyp_select<-"H0"
source("scripts/oncology_batches_1-1.R")

#Large n H1
scen_select<-"large_n"
hyp_select<-"H1"
source("scripts/oncology_batches_1-1.R")

#Large n H0
scen_select<-"large_n"
hyp_select<-"H0"
source("scripts/oncology_batches_1-1.R")

#Extra IPCW H1 (this only uses small n)
scen_select<-"IPCW_extra"
hyp_select<-"H1"
source("scripts/oncology_batches_1-1.R")

#Extra IPCW H0 (this only uses small n)
scen_select<-"IPCW_extra"
hyp_select<-"H0"
source("scripts/oncology_batches_1-1.R")

