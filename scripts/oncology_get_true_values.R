#Get true values

# devtools::install()
# renv::restore()
# devtools::document()
rm(list=ls())

devtools::load_all()
library(CI.RCT.Sim)
library(parallel)
library(survival)

#Make scenario table, the code now also adds the column with scneario short names
save_param_tab<-FALSE #need to create folder "data" if not present, as it is saved there
source("scripts/oncology_make_scenario_table_1-1.R")

#Run different scenario settings:
#need to create folder "results" if not present!!
if(!file.exists("results"))  dir.create("results")

#Set number of iterations

pre_N_sim <- 20
ev_soll_for_true_value<-10000


#SimClean()

#scen_tab <- readxl::read_xlsx("data/oncology_scenario_list.xlsx")
n_scenarios<-dim(scen_tab)[1]

scen_set<-(1:n_scenarios)[grepl("H1",scen_tab$block_nam)]
result_name_note<-"true_values"


# Alpha
alpha <- 0.05 #two sided, tests will be one-sided using alpha/2, confidence intervals are mostly hard coded to 0.95


###################################

# Derive true treatment effect -------------------------------------------

sim_parameters <- oncology_scenario() |>
  oncology_scenario_set_truevalues()




#pre_sim_parameters <- oncology_scenario()


pre_my_analyse <- list(
  truth = function(condition, dat, fixed_objects = NULL) {
    mod <- survival::coxph(Surv(time = event_time, event = ev) ~ trt + X_0 + W_0, data = dat)
    HR <- exp(coef(mod)[1])
    list(HR = HR)
  }
)

pre_my_summarise <- create_summarise_function(
  truth = summarise_estimator( #test if this can be an object defined before
    est = HR,
    real = exp(beta_death[[1]][6]), #actually not required, but we need a placeholder here
    null = 1
  )
)

summy<-summarise_estimator( #test if this can be an object defined before
  est = HR,
  real = exp(beta_death[[1]][6]), #actually not required, but we need a placeholder here
  null = 1
)
pre_my_summarise <- create_summarise_function(
  truth = summy
)

cl <- makeCluster(detectCores(logical = FALSE) - 1)
clusterEvalQ(cl, {
  library("CI.RCT.Sim")
  #library("survival")
})

#SimClean()

pre_results <- runSimulation(
  design = sim_parameters[scen_set,],
  replications = pre_N_sim,
  generate = generate_oncology,
  analyse = pre_my_analyse,
  summarise = pre_my_summarise,
  fixed_objects = list(allow_switch = FALSE, logHR_assumed = NULL, ev_soll = ev_soll_for_true_value, allow_random_cens = TRUE),
  parallel = TRUE,
  cl = cl
)
as.data.frame(pre_results)
stopCluster(cl)
pre_results$truth.mean_est
names(pre_results)
out<-data.frame(scen=scen_set,trueHR=pre_results$truth.mean_est)
out
filename<-paste("scripts/trueHR_nsim",pre_N_sim,"_n_ev",ev_soll_for_true_value,".csv",sep="")
write.csv(out,file=filename)

# Save results -----------------------------------------------------------
#path="results/"
#file_name<-paste(path,result_name_note,"_",hyp_select,"_",format(Sys.time(), paste0("results_onco_","nsim",N_sim,"_", Sys.info()["nodename"], "%Y-%m-%d_%H%M.Rdata")),sep="")
#file_name
#save(results, main_sessioninfo, nodes_sessioninfo, file = file_name)


#results
#A<-as.data.frame(results)
#rej<-grepl("test.rejection_0.025",names(A))
#cover<-grepl("est.coverage",names(A))
#A[rej]
#A[cover]
