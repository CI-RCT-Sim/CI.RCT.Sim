# Variant of test_ps_coverage.R that swaps in a robust (sandwich) vcov.
# We pass vcov. = sandwich::vcovHC(., type = "HC0") to emmeans, which causes
# the entire pipeline (regrid, contrast, confint, test) to use that vcov.
#
# Note: this still does NOT account for uncertainty in the estimated
# principal score (weights treated as fixed). It only fixes the part where
# glm's model-based SE is wrong for non-frequency weights.

suppressPackageStartupMessages({
  library(dplyr)
  library(emmeans)
  library(sandwich)
})

set.seed(1)

simulate_one <- function(n = 4000,
                         p_C_trt1 = 0.7,
                         p_C_trt0 = 0.7,
                         risk_ctrl_comp = 0.05,
                         VE_true = 0.3) {
  n_per <- n / 2
  trt <- c(rep(1, n_per), rep(0, n_per))
  C  <- rbinom(n, 1, ifelse(trt == 1, p_C_trt1, p_C_trt0))
  risk <- ifelse(trt == 1 & C == 1, risk_ctrl_comp * (1 - VE_true),
                 ifelse(trt == 0 & C == 1, risk_ctrl_comp,
                        0.05))
  evt <- rbinom(n, 1, risk)
  V <- rbinom(n, 1, 0.3)
  W <- rbinom(n, 1, 0.3)
  data.frame(trt = trt, C = C, evt = evt, V = V, W = W)
}

# Fixed analysis using sandwich vcov ------------------------------------
analyse_sandwich <- function(dat, ci_level = 0.95, VE_margin = 0.3) {
  dat1 <- dat |>
    within({
      C   <- factor(C,   levels = c("0", "1"))
      trt <- factor(trt, levels = c("1", "0"))
    })
  mod_ps <- glm(C ~ V + W, subset = (trt == 1), data = dat1, family = binomial())
  odds <- exp(predict(mod_ps, newdata = dat1, type = "link"))

  dat1 <- dat1 |>
    mutate(
      weight = case_when(
        (trt == 1) & (C == "1") ~ 1,
        (trt == 0) & (C == "0") ~ 0,
        (trt == 1) & (C == "0") ~ 0,
        (trt == 0) & (C == "1") ~ odds
      )
    )

  outcome_mod <- suppressWarnings(
    glm(evt ~ trt + V + W, weights = weight, family = binomial(), data = dat1)
  )

  # Key change: pass robust vcov to emmeans. The function vcovHC from sandwich
  # gives the heteroskedasticity-consistent (Huber/White) sandwich estimator,
  # which is the right SE for a weighted GLM treating weights as known.
  emm  <- emmeans(outcome_mod, ~ trt, vcov. = function(x) vcovHC(x, type = "HC0"))
  lemm <- regrid(emm, "log")
  ctr  <- contrast(lemm, method = "pairwise", type = "response")

  ci_two <- confint(ctr, level = ci_level)             # two-sided 95% CI
  tst    <- test(ctr, null = log(1 - VE_margin), side = "<")  # one-sided test

  list(
    p        = tst$p.value,
    VE       = 1 - ci_two$ratio,
    VE_lower = 1 - ci_two$asymp.UCL,
    VE_upper = 1 - ci_two$asymp.LCL
  )
}

# Naive (model-based) version for comparison ----------------------------
analyse_naive <- function(dat, ci_level = 0.95, VE_margin = 0.3) {
  dat1 <- dat |>
    within({
      C   <- factor(C,   levels = c("0", "1"))
      trt <- factor(trt, levels = c("1", "0"))
    })
  mod_ps <- glm(C ~ V + W, subset = (trt == 1), data = dat1, family = binomial())
  odds <- exp(predict(mod_ps, newdata = dat1, type = "link"))

  dat1 <- dat1 |>
    mutate(
      weight = case_when(
        (trt == 1) & (C == "1") ~ 1,
        (trt == 0) & (C == "0") ~ 0,
        (trt == 1) & (C == "0") ~ 0,
        (trt == 0) & (C == "1") ~ odds
      )
    )

  outcome_mod <- suppressWarnings(
    glm(evt ~ trt + V + W, weights = weight, family = binomial(), data = dat1)
  )

  emm  <- emmeans(outcome_mod, ~ trt)
  lemm <- regrid(emm, "log")
  ctr  <- contrast(lemm, method = "pairwise", type = "response")
  ci_two <- confint(ctr, level = ci_level)
  tst    <- test(ctr, null = log(1 - VE_margin), side = "<")

  list(
    p        = tst$p.value,
    VE       = 1 - ci_two$ratio,
    VE_lower = 1 - ci_two$asymp.UCL,
    VE_upper = 1 - ci_two$asymp.LCL
  )
}

# Monte Carlo -----------------------------------------------------------
B <- 2000
VE_true <- 0.3

cat("Running", B, "replicates under H0 (VE_AC = 0.3) ...\n")

res_naive    <- matrix(NA_real_, B, 3, dimnames = list(NULL, c("p","lo","hi")))
res_sandwich <- res_naive

pb <- txtProgressBar(min = 0, max = B, style = 3)
for (b in seq_len(B)) {
  d  <- simulate_one(VE_true = VE_true)
  rn <- tryCatch(analyse_naive(d),    error = function(e) NULL)
  rs <- tryCatch(analyse_sandwich(d), error = function(e) NULL)
  if (!is.null(rn)) res_naive[b, ]    <- c(rn$p, rn$VE_lower, rn$VE_upper)
  if (!is.null(rs)) res_sandwich[b, ] <- c(rs$p, rs$VE_lower, rs$VE_upper)
  setTxtProgressBar(pb, b)
}
close(pb)

summarise_results <- function(M, label, VE_true) {
  ok <- complete.cases(M)
  p  <- M[ok, "p"]; lo <- M[ok, "lo"]; hi <- M[ok, "hi"]
  cat(sprintf("\n--- %s (n_eval = %d / %d) ---\n", label, sum(ok), nrow(M)))
  cat(sprintf("  Rejection at one-sided alpha = 0.025: %.4f  (target 0.025)\n", mean(p < 0.025)))
  cat(sprintf("  Rejection at one-sided alpha = 0.05 : %.4f  (target 0.05)\n",  mean(p < 0.05)))
  cat(sprintf("  CI coverage of VE_true=%.2f         : %.4f  (target 0.95)\n", VE_true,
              mean((lo <= VE_true) & (VE_true <= hi))))
  cat(sprintf("    miss low  (VE_true < VE_lower)   : %.4f\n", mean(VE_true < lo)))
  cat(sprintf("    miss high (VE_true > VE_upper)   : %.4f\n", mean(VE_true > hi)))
  cat(sprintf("  Mean VE_upper                       : %.4f\n", mean(hi)))
  cat(sprintf("  Mean VE_lower                       : %.4f\n", mean(lo)))
}

summarise_results(res_naive,    "NAIVE model-based vcov", VE_true)
summarise_results(res_sandwich, "SANDWICH vcovHC(type='HC0')", VE_true)
