# List of analysis functions ---------------------------------------------

my_analyse <- list(
  ## Treatment policy estimands
  ipwtp = analyse_diabetes_ipw(strategy = "treatment_policy"),
  mmrmtp = analyse_diabetes_mmrm(strategy = "treatment_policy"),
  mitp = analyse_diabetes_mi(strategy = "treatment_policy"),
  ## Hypothetical estimands
  ipwhyp = analyse_diabetes_ipw(strategy = "hypothetical"),
  dmhyp = analyse_diabetes_demediation(),
  gcomhyp = analyse_diabetes_gcomputation(),
  mmrmhyp = analyse_diabetes_mmrm(strategy = "hypothetical"),
  mihyp = analyse_diabetes_mi(strategy = "hypothetical"),
  describe = function(condition, dat, fixed_objects = NULL) {
    result <- list(
      n_pat  = nrow(dat),
      n_ctrl = sum(dat$trt == 0),
      n_trt  = sum(dat$trt == 1),
      n_resc = sum(dat$rescue_start < 14, na.rm = TRUE),
      n_miss = sum(dat$m_start < 13, na.rm = TRUE),
      naive  = mean(dat[which(dat$trt == 1),]$y12-dat[which(dat$trt == 1),]$y0, na.rm = TRUE) -
        mean(dat[which(dat$trt == 0),]$y12-dat[which(dat$trt == 0),]$y0, na.rm = TRUE)
    )
    result
  }

)

my_analyse <- wrap_all_in_trycatch(my_analyse)

# List of summarisation functions ----------------------------------------
# summarise_estimator and summarise_test are generic summarisation
# functions from CI.RCT.Sim / SimDesign

sum_tp <- summarise_estimator(est = coef, real = tp_eff, lower = ci_lower,
  upper = ci_upper, null = 0, name = "est")

sum_hyp <- summarise_estimator(est = coef, real = eff_true, lower = ci_lower,
  upper = ci_upper, null = 0, name = "est")

sum_test <- summarise_test(alpha, name = "test")


my_summarise <- create_summarise_function(
  # bias, SD, coverage etc. for the treatment effect at final visit
  ## Treatment policy estimands
  ipwtp = sum_tp,
  mmrmtp = sum_tp,
  mitp = sum_tp,
  ## Hypothetical estimands
  ipwhyp = sum_hyp,
  dmhyp = sum_hyp,
  gcomhyp = sum_hyp,
  mmrmhyp = sum_hyp,
  mihyp = sum_hyp,
  # rejection rates
  ## Treatment policy estimands
  ipwtp = sum_test,
  mmrmtp = sum_test,
  mitp = sum_test,
  ## Hypothetical testimands
  ipwhyp = sum_test,
  dmhyp = sum_test,
  gcomhyp = sum_test,
  mmrmhyp = sum_test,
  mihyp = sum_test,
  mmrmtp = function(condition, results, fixed_objects = NULL){
    data.frame(n_conv = sum(results$converged),
               fallbacks = sum(results$fallback),
               unstruc = sum(results$covariance == "us"))},
  mmrmhyp = function(condition, results, fixed_objects = NULL){
    data.frame(n_conv = sum(results$converged),
               fallbacks = sum(results$fallback),
               unstruc = sum(results$covariance == "us"))},
  describe = summarise_describe()
)
