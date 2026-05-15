
#to be run from oncology_run_all.R

# devtools::install()
# renv::restore()
# devtools::document()
#devtools::load_all()
#rm(list=ls())
#library(CI.RCT.Sim)
#library(parallel)
#library(survival)

#SimClean()

#scen_tab <- readxl::read_xlsx("data/oncology_scenario_list.xlsx")
n_scenarios<-dim(scen_tab)[1]


# Settings to calculate true value
#pre_N_sim <- 2#10
#ev_soll_for_true_value<-100#0

# Iterations and scenarios
#N_sim <- 5
#scen_set<-c(52,57,64)#53 #H0,, #1:3
##scen_select<-"all"

#### run with these three settings separately:
#scen_select<-"small_n"
#scen_select<-"large_n"
# scen_select<-"IPCW_extra"

#hyp_select<-"H1"
#hyp_select<-"H0"
#

set_hyp<-grepl(hyp_select,scen_tab$block_nam)

result_name_note<-""

if(scen_select=="all") {
  scen_set<-1:n_scenarios
  result_name_note<-"scen_all"
}
if(scen_select=="small_n") { #small n is high effect size
  scen_set<-(1:n_scenarios)[grepl("high",scen_tab$block_nam) & set_hyp]
  result_name_note<-"scen_small_n"
}
if(scen_select=="large_n") {
  scen_set<-(1:n_scenarios)[grepl("low",scen_tab$block_nam) & set_hyp]
  result_name_note<-"scen_large_n"
}
if(scen_select=="IPCW_extra") {
  scen_set<-(1:n_scenarios)[grepl("high",scen_tab$block_nam) & set_hyp & (grepl("Core",scen_tab$scenario_name) | grepl("random censoring",scen_tab$scenario_name))]
  result_name_note<-"scen_IPCW_extra"
}



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
#stopCluster(cl)

# Under H0 the true effect is HR = 1
pre_results[which(sapply(pre_results$beta_death, `[[`, 6) == 0),]$truth.mean_est <- 1

# Define parameter values and derived quantities -------------------------

#sim_parameters <- oncology_scenario() |>
#  oncology_scenario_set_truevalues() |>
#  dplyr::mutate(true_eff = pre_results$truth.mean_est) #redundant, pool in first step

#sim_parameters <- sim_parameters |>
#  dplyr::mutate(true_eff = pre_results$truth.mean_est)
sim_parameters$true_eff<-NA
sim_parameters$true_eff[scen_set]<-pre_results$truth.mean_est
# Constants for simulation -----------------------------------------------



# List of analysis functions ---------------------------------------------

#all or large n case
analysis_functions_list<-list(
  rpsftm_rc = analyse_oncology_rpsftm(recensor = TRUE),
  rpsftm = analyse_oncology_rpsftm(recensor = FALSE),
  tse_rc = analyse_oncology_TSE(recensor = TRUE),
  tse = analyse_oncology_TSE(recensor = FALSE),
  #gformula = analyse_oncology_gformula(B = 20),
  ipw = analyse_oncology_ipw(),
  #itt = analyse_oncology_itt(),
  cens = analyse_oncology_cens()
)

if(scen_select=="small_n") analysis_functions_list<-c(analysis_functions_list,gformula = analyse_oncology_gformula(B = 20))


if(scen_select=="IPCW_extra") {
  #IPCW_functions
  analysis_functions_list<-list(
    rpsftm = analyse_oncology_rpsftm(recensor = FALSE),
    tse = analyse_oncology_TSE(recensor = FALSE),
    gformula = analyse_oncology_gformula(B = 20),
    ipw = analyse_oncology_ipw(),

    rpsftm_IPCW =  analyse_oncology_mixed(method="RPSFTM",recensor = TRUE,B = 100,trunc_weights = 5,use_censoring_IPW = TRUE,requ_n_cens = 5),
    tse_IPCW =  analyse_oncology_mixed(method="TSE",recensor = TRUE,B = 100,trunc_weights = 5,use_censoring_IPW = TRUE,requ_n_cens = 5),
    gformula_IPCW = analyse_oncology_gformula(B = 20,use_censoring_IPW=TRUE, requ_n_cens=5, trunc_weights=5),
    ipw_IPCW =analyse_oncology_ipw2(use_censoring_IPW = TRUE, trunc_weights = 5, requ_n_cens = 5)
  )

}

my_analyse <- c(
  analysis_functions_list,
  list(
  #rpsftm_rc = analyse_oncology_rpsftm(recensor = TRUE),
  #rpsftm = analyse_oncology_rpsftm(recensor = FALSE),
  #tse_rc = analyse_oncology_TSE(recensor = TRUE),
  #tse = analyse_oncology_TSE(recensor = FALSE),
  #gformula = analyse_oncology_gformula(B = 20),
  #ipw = analyse_oncology_ipw(),
  #itt = analyse_oncology_itt(),
  #cens = analyse_oncology_cens(),
  describe = function(condition, dat, fixed_objects = NULL) {
    tabulate_helper <- function(dat, var) {
      tmp <- list(
        sum(dat[, var]),
        sum(dat[dat$trt == 0, var]),
        sum(dat[dat$trt == 1, var])
      )

      names(tmp) <- c(var, paste0(var, "_ctrl"), paste0(var, "_trt"))
      tmp
    }
    result <- list(
      n_pat = nrow(dat),
      n_ctrl = sum(dat$trt == 0),
      n_trt = sum(dat$trt == 1),
      n_switch = sum(dat$switch),
      max_followup = max(dat$event_time),
      sufficient_events = sum(dat$ev) >= condition$ev_soll,
      n_random_cens = sum(dat$random_cens),
      sufficient_random_cens = sum(dat$random_cens)>=5
    )
    result <- c(result, tabulate_helper(dat, "ev"))
    if (!is.null(attr(dat, "followup"))) {
      result$study_time <- attr(dat, "followup")
    } else {
      result$study_time <- NA_real_
    }
    if (hasName(dat, "ice")) {
      result <- c(result, tabulate_helper(dat, "ice"))
    }
    if (hasName(dat, "subgroup")) {
      result <- c(result, tabulate_helper(dat, "subgroup"))
    }
    result
  }
))

my_analyse <- wrap_all_in_trycatch(my_analyse)

# List of summarisation functions ----------------------------------------
# summarise_estimator and summarise_test are generic summarisation
# functions from CI.RCT.Sim / SimDesign


summy<-summarise_estimator(
  est = HR,
  real = true_eff,
  lower = low,
  upper = up,
  null = 1,
  name = "est"
)

summarise_test_one_sided_for_HR<-function (alpha, name = NULL) {
  res <- function(condition, results, fixed_objects) {
    results$p<-results$p/2
    results$p[results$HR>1]<-1-results$p[results$HR>1]
    rejection_tmp <- setNames(as.data.frame(as.list(colMeans(outer(results$p,
                                                                   alpha, FUN = `<`), na.rm = TRUE))), paste0("rejection_",
                                                                                                              alpha))
    missing_tmp <- setNames(as.data.frame(as.list(colSums(outer(results$p,
                                                                1 - alpha, FUN = function(p, a) {
                                                                  is.na(p)
                                                                })))), paste0("N_missing_", alpha))
    results_tmp <- cbind(rejection_tmp, missing_tmp, N = nrow(results))
    results_tmp$mean_n_pat <- NA_real_
    results_tmp$sd_n_pat <- NA_real_
    results_tmp$mean_n_evt <- NA_real_
    results_tmp$sd_n_evt <- NA_real_
    results_tmp$N_missing_n_pat <- NA_real_
    results_tmp$N_missing_n_evt <- NA_real_
    if (hasName(results, "N_pat")) {
      results_tmp$mean_n_pat <- mean(results$N_pat, na.rm = TRUE)
      results_tmp$sd_n_pat <- sd(results$N_pat, na.rm = TRUE)
      results_tmp$N_missing_n_pat <- sum(is.na(results$N_pat))
    }
    if (hasName(results, "N_evt")) {
      results_tmp$mean_n_evt <- mean(results$N_evt, na.rm = TRUE)
      results_tmp$sd_n_evt <- sd(results$N_evt, na.rm = TRUE)
      results_tmp$N_missing_n_evt <- sum(is.na(results$N_evt))
    }
    results_tmp
  }
  attr(res, "name") <- name
  res
}

sumtest<-summarise_test_one_sided_for_HR(
  alpha/2,
  name = "test"
)

my_summarise <- create_summarise_function(
  # bias, SD, coverage etc. for the treatment effect at final visit
  rpsftm_rc = summy,
  rpsftm = summy,
  tse_rc = summy,
  tse = summy,
  gformula = summy,
  ipw = summy,
  itt = summy,
  cens = summy,
  # rejection rates
  rpsftm_rc = sumtest,
  rpsftm = sumtest,
  tse_rc = sumtest,
  tse = sumtest,
  gformula = sumtest,
  ipw = sumtest,
  itt = sumtest,
  cens = sumtest,
  describe = summarise_describe()
)

# Run the simulations ----------------------------------------------------

#cl <- makeCluster(detectCores(logical = FALSE) - 1)
#clusterEvalQ(cl, {
#  library("CI.RCT.Sim")
#})
clusterExport(cl = cl, varlist = c("alpha"))

main_sessioninfo <- sessionInfo()
nodes_sessioninfo <- clusterEvalQ(cl, {
  sessionInfo()
})

results <- runSimulation(
  design = sim_parameters[scen_set,],
  replications = N_sim,
  generate = generate_oncology,
  analyse = my_analyse,
  summarise = my_summarise,
  fixed_objects = list(allow_switch = TRUE, logHR_assumed = NULL, ev_soll = NULL, allow_random_cens = TRUE),
  parallel = TRUE,
  cl = cl
)

stopCluster(cl)

# Save results -----------------------------------------------------------
path="results/"
file_name<-paste(path,result_name_note,"_",hyp_select,"_",format(Sys.time(), paste0("results_onco_","nsim",N_sim,"_", Sys.info()["nodename"], "%Y-%m-%d_%H%M.Rdata")),sep="")
file_name
save(results, main_sessioninfo, nodes_sessioninfo, file = file_name)


#results
A<-as.data.frame(results)
rej<-grepl("test.rejection_0.025",names(A))
cover<-grepl("est.coverage",names(A))
#A[rej]
#A[cover]
