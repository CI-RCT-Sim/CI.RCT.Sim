# List of analysis functions ---------------------------------------------

my_analyse <- list(
  # both V and W observed
  iv       = analyse_vaccine_ivreg(ci_level = 1-alpha_ci, VE_margin = 0.3),
  iv2      = analyse_vaccine_ivreg2(ci_level = 1-alpha_ci, VE_margin = 0.3),
  ps_cov   = analyse_vaccine_ps(ci_level = 1-alpha_ci, VE_margin = 0.3, covariates_in_outcomes_model = TRUE),
  ps_nocov = analyse_vaccine_ps(ci_level = 1-alpha_ci, VE_margin = 0.3, covariates_in_outcomes_model = FALSE),
  pp       = analyse_vaccine_pp(ci_level = 1-alpha_ci, VE_margin = 0.3),
  # with trt x W interaction
  ps_cov_winter   = analyse_vaccine_ps(ci_level = 1-alpha_ci, VE_margin = 0.3, covariates_in_outcomes_model = TRUE, W_interaction = TRUE),
  ps_nocov_winter = analyse_vaccine_ps(ci_level = 1-alpha_ci, VE_margin = 0.3, covariates_in_outcomes_model = FALSE, W_interaction = TRUE),
  pp_winter       = analyse_vaccine_pp(ci_level = 1-alpha_ci, VE_margin = 0.3, W_interaction = TRUE),
  # V unobserved
  iv_vunobs       = analyse_vaccine_ivreg(ci_level = 1-alpha_ci, VE_margin = 0.3, V_unobserved=TRUE),
  iv_vunobs2      = analyse_vaccine_ivreg2(ci_level = 1-alpha_ci, VE_margin = 0.3, V_unobserved=TRUE),
  ps_cov_vunobs   = analyse_vaccine_ps(ci_level = 1-alpha_ci, VE_margin = 0.3, covariates_in_outcomes_model = TRUE, V_unobserved=TRUE),
  ps_nocov_vunobs = analyse_vaccine_ps(ci_level = 1-alpha_ci, VE_margin = 0.3, covariates_in_outcomes_model = FALSE, V_unobserved=TRUE),
  pp_vunobs       = analyse_vaccine_pp(ci_level = 1-alpha_ci, VE_margin = 0.3, V_unobserved=TRUE),
  # with trt x W interaction
  ps_cov_vunobs_winter   = analyse_vaccine_ps(ci_level = 1-alpha_ci, VE_margin = 0.3, covariates_in_outcomes_model = TRUE, V_unobserved=TRUE, W_interaction = TRUE),
  ps_nocov_vunobs_winter = analyse_vaccine_ps(ci_level = 1-alpha_ci, VE_margin = 0.3, covariates_in_outcomes_model = FALSE, V_unobserved=TRUE, W_interaction = TRUE),
  pp_vunobs_winter       = analyse_vaccine_pp(ci_level = 1-alpha_ci, VE_margin = 0.3, V_unobserved=TRUE, W_interaction = TRUE),
  # W unobserved
  iv_wunobs       = analyse_vaccine_ivreg(ci_level = 1-alpha_ci, VE_margin = 0.3, W_unobserved=TRUE),
  iv_wunobs2      = analyse_vaccine_ivreg2(ci_level = 1-alpha_ci, VE_margin = 0.3, W_unobserved=TRUE),
  ps_cov_wunobs   = analyse_vaccine_ps(ci_level = 1-alpha_ci, VE_margin = 0.3, covariates_in_outcomes_model = TRUE, W_unobserved=TRUE),
  ps_nocov_wunobs = analyse_vaccine_ps(ci_level = 1-alpha_ci, VE_margin = 0.3, covariates_in_outcomes_model = FALSE, W_unobserved=TRUE),
  pp_wunobs       = analyse_vaccine_pp(ci_level = 1-alpha_ci, VE_margin = 0.3, W_unobserved=TRUE),
  # both V and W unobserved
  iv_vwunobs       = analyse_vaccine_ivreg(ci_level = 1-alpha_ci, VE_margin = 0.3, V_unobserved=TRUE, W_unobserved=TRUE),
  iv_vwunobs2      = analyse_vaccine_ivreg2(ci_level = 1-alpha_ci, VE_margin = 0.3, V_unobserved=TRUE, W_unobserved=TRUE),
  ps_cov_vwunobs   = analyse_vaccine_ps(ci_level = 1-alpha_ci, VE_margin = 0.3, covariates_in_outcomes_model = TRUE, V_unobserved=TRUE, W_unobserved=TRUE),
  ps_nocov_vwunobs = analyse_vaccine_ps(ci_level = 1-alpha_ci, VE_margin = 0.3, covariates_in_outcomes_model = FALSE, V_unobserved=TRUE, W_unobserved=TRUE),
  pp_vwunobs       = analyse_vaccine_pp(ci_level = 1-alpha_ci, VE_margin = 0.3, V_unobserved=TRUE, W_unobserved=TRUE)
)

my_analyse <- wrap_all_in_trycatch(my_analyse)

# List of summarisation functions ----------------------------------------
# summarise_estimator and summarise_test are generic summarisation
# functions from CI.RCT.Sim / SimDesign

# those functions should be used for summarisation for all the methods
tmp_functions <- list(
  summarise_estimator(VE, VE, VE_lower, VE_upper, null=0.3, name="est"),
  summarise_test(alpha_test, name="test"),
  summarise_estimator(VE_sandwich, VE, VE_lower_sandwich, VE_upper_sandwich, null=0.3, name="est_sandwich")
)

# create list with summarisation functions as elements and names of methods as names
# call create_summarise_function with thoses
my_summarise <- lapply(tmp_functions, \(fn){
  lapply(names(my_analyse), \(n){
    fn
  }) |> setNames(names(my_analyse))
}) |>
  do.call(c, args=_) |>
  do.call(create_summarise_function, args=_)

