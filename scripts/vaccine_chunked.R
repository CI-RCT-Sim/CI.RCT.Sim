# devtools::install()
# renv::restore()
library(CI.RCT.Sim)
library(parallel)

source("scripts/vaccine_scenario_classes.R")
source("scripts/vaccine_analyse_summarise.R")

# Define parameter values and derived quantities -------------------------

scenario <- Sys.getenv("scenario")

message(paste("Scenario", scenario, "selected"))

selected_scenario <- switch(
  scenario,
  A1 = vaccine_scenario_A1,
  A2 = vaccine_scenario_A2,
  B1 = vaccine_scenario_B1,
  C1 = vaccine_scenario_C1,
  D1 = vaccine_scenario_D1,
  extra = vaccine_scenario_extra,
  stop("unknown scenario selected or selection missing")
)


# setup chunks of scenarios -----------------------------------------------

start_chunk <- Sys.getenv("start") |>
  strtoi()

chunksize <- 5
N_scenarios <- selected_scenario |>
  nrow()
N_chunks <- ceiling(N_scenarios / chunksize)

if(start_chunk > N_chunks){
  stop(paste0("Start chunk, ", start_chunk, " is larger than number of chunks ", N_chunks))
}

rows_chunks <- lapply(1:N_chunks, \(i){
  rows <- 1:chunksize + (i-1)*chunksize
  rows <- rows[rows <= N_scenarios]
  rows
})

# run simulations ---------------------------------------------------------

for(i in start_chunk:length(rows_chunks)){
  message(paste0("Running chunk ", i, " of ", N_chunks, ", rows ", paste0(rows_chunks[[i]], collapse=", ")))

  sim_parameters <- selected_scenario |>
    vaccine_scenario_set_beta_A1_relative() |>
    vaccine_scenario_set_gamma_0() |>
    vaccine_scenario_set_true_eff() |>
    vaccine_scenario_set_samplesize() |>
    within({
      VE = 1-rr_ps
      scenario_nr = seq_along(VE)
    }) |>
    _[rows_chunks[[i]], ]

  message(paste(nrow(sim_parameters), "rows"))

  # Constants for simulation -----------------------------------------------

  N_sim <- 10000
  alpha_ci <- 0.05
  alpha_test <- c(0.05, 0.025)

  # Run the simulations ----------------------------------------------------

  message(paste(length(my_analyse), "analysis functions"))
  message(paste(length(environment(my_summarise)$summarise_functions), "summarise functions"))
  message(paste("setting up cluster with", detectCores(logical=TRUE)-1 , "cores"))

  cl <- makeCluster(detectCores(logical=TRUE)-1)
  clusterEvalQ(cl, {
    library("CI.RCT.Sim")
  })

  clusterExport(cl = cl, varlist = c("alpha_ci", "alpha_test"))

  main_sessioninfo <- sessionInfo()
  nodes_sessioninfo <- clusterEvalQ(cl, {
    sessionInfo()
  })

  message(paste("running simulations,", nrow(sim_parameters), "scenarios,", N_sim, "replications"))

  results <- runSimulation(
    design = sim_parameters,
    replications = N_sim,
    generate = generate_vaccine,
    analyse = my_analyse,
    summarise = my_summarise,
    fixed_objects = list(include_unobserved=FALSE),
    parallel = TRUE,
    cl = cl
  )

  message("stopping cluster")
  stopCluster(cl)

  # Save results -----------------------------------------------------------

  out_filename <- format(Sys.time(), paste0("results_vaccine_scenario_", scenario, "_chunk_", i, "_", Sys.info()["nodename"], "%Y-%m-%d_%H%M.Rdata"))
  message(paste("saving results to", out_filename))
  save(results, main_sessioninfo, nodes_sessioninfo, file=out_filename)

}

