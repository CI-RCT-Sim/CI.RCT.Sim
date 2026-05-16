
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


n_scenarios<-dim(scen_tab)[1]

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


#make batches of up to batch_length scenarios
#scen_set
n_selected_scen<-length(scen_set)
n_batches<-ceiling(n_selected_scen/batch_length)
n_batches
size_last_batch<-n_selected_scen%%batch_length
if(n_batches==1) {
  batch_list<-list(scen_set)
} else {

  batch_list<-vector(length=n_batches,mode="list")


  for(vv in 1:(n_batches-1)) batch_list[[vv]]<-scen_set[1:batch_length+batch_length*(vv-1)]
  batch_list[[length(batch_list)]]<-scen_set[batch_length*(n_batches-1)+1:size_last_batch]

}
#batch_index<-1
for(batch_index in 1:n_batches) {
  scen_set<-batch_list[[batch_index]]
  batch_name<-paste("sc","_",scen_set[1],"to",scen_set[length(scen_set)],sep="")

  sim_parameters <- oncology_scenario() |>
    oncology_scenario_set_truevalues()





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

cl <- makeCluster(detectCores(logical = FALSE) - 1)
clusterEvalQ(cl, {
  library("CI.RCT.Sim")
  #library("survival")
})
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
file_name<-paste(path,result_name_note,"_",hyp_select,"_",batch_name,"_",format(Sys.time(), paste0("results_onco_","nsim",N_sim,"_", Sys.info()["nodename"], "%Y-%m-%d_%H%M.Rdata")),sep="")
file_name
save(results, main_sessioninfo, nodes_sessioninfo, file = file_name)

}

#results
#A<-as.data.frame(results)
#rej<-grepl("test.rejection_0.025",names(A))
#cover<-grepl("est.coverage",names(A))
#A[rej]
#A[cover]
