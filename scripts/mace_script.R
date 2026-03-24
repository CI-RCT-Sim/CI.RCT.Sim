devtools::load_all()

# -------------------------------------------------------------------
# Define parameter values and derived quantities
# -------------------------------------------------------------------

sim_parameters <- mace_scenario() |> true_trt_mace()

# -------------------------------------------------------------------
# Constants for simulation
# -------------------------------------------------------------------

N_sim <- 100

# -------------------------------------------------------------------
# List of analysis functions
# -------------------------------------------------------------------

my_analyse <- list(
  cox_nocov=analyse_mace_cox_nocov,
  cox_cov=analyse_mace_cox_cov,
  ipw_nocov=analyse_mace_ipw_nocov,
  ipw_cov=analyse_mace_ipw_cov
)


my_summarise <- create_summarise_function(
  cox_nocov=summarise_estimator(est=coef, real=true_trt, est_sd=se, name="est"),
  cox_cov  =summarise_estimator(est=coef, real=true_trt, est_sd=se, name="est"),
  ipw_nocov=summarise_estimator(est=coef, real=true_trt, est_sd=se, name="est"),
  ipw_cov  =summarise_estimator(est=coef, real=true_trt, est_sd=se, name="est"),
  cox_nocov=summarise_estimator(est=HR, real=exp(true_trt), lower=CI.l, upper=CI.u, null=1, name="HR"),
  cox_cov  =summarise_estimator(est=HR, real=exp(true_trt), lower=CI.l, upper=CI.u, null=1, name="HR"),
  ipw_nocov=summarise_estimator(est=HR, real=exp(true_trt), lower=CI.l, upper=CI.u, null=1, name="HR"),
  ipw_cov  =summarise_estimator(est=HR, real=exp(true_trt), lower=CI.l, upper=CI.u, null=1, name="HR"),
  cox_nocov=summarise_test(alpha=0.05, name="test"),
  cox_cov  =summarise_test(alpha=0.05, name="test"),
  ipw_nocov=summarise_test(alpha=0.05, name="test"),
  ipw_cov  =summarise_test(alpha=0.05, name="test")
)

# -------------------------------------------------------------------
# Run the simulations
# -------------------------------------------------------------------

results <- runSimulation(
  design = sim_parameters,
  replications = N_sim,
  generate = generate_mace,
  analyse = my_analyse,
  summarise = my_summarise
)

# -------------------------------------------------------------------
# Inspect results
# -------------------------------------------------------------------

methods <- c("cox_nocov", "cox_cov", "ipw_nocov", "ipw_cov")

# estimators HR scale
results[, paste(methods, "HR", "bias", sep=".")]
results[, paste(methods, "HR", "width", sep=".")]

# estimator log-HR scale
results[, paste(methods, "est", "bias", sep=".")]
results[, paste(methods, "est", "sd_est", sep=".")]

# test decision from CI and test
results[, paste(methods, "test", "rejection_0.05", sep=".")]
1-results[, paste(methods, "HR", "null_cover", sep=".")]
