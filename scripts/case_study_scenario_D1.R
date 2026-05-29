set.seed(123)

source("scripts/vaccine_scenario_classes.R")

condition <- vaccine_scenario_D1 |>
  _[1,] |>
  vaccine_scenario_set_beta_A1_relative() |>
  vaccine_scenario_set_gamma_0() |>
  vaccine_scenario_set_true_eff() |>
  vaccine_scenario_set_samplesize()

dat <- generate_vaccine(condition)

1-exp(condition$beta_A2)
1-condition$rr_ps

# res <- list(
#   analyse_vaccine_ps(covariates_in_outcomes_model = TRUE)(condition, dat),
#   analyse_vaccine_ps(covariates_in_outcomes_model = FALSE)(condition, dat),
#   analyse_vaccine_ivreg()(condition, dat),
#   analyse_vaccine_pp()(condition, dat)
# ) |>
#   lapply(as.data.frame) |>
#   purrr::list_rbind() |>
#   _[c("VE", "VE_lower", "p")]
#
# (res$VE * 100) |>
#   round(1) |>
#   cat(sep="\n")

# IV regression -----------------------------------------------------------
library(emmeans)

# prepare data
dat1  <- dat |>
  within({
    # Indicator 1 if complier in treatment group
    T <- factor(ifelse(((trt==1) & (C==1)), 1L, 0L), levels=c("1", "0"))
    # recode 0/1 coded variables as factors
    V <- factor(V)
    W <- factor(W)
    trt <- factor(trt)
    # recode F/T coded variable to 0/1
    evt <- as.integer(evt)
  })

# estimate stage 1 and save to dat1
x_stage1 <- model.matrix(~ trt + V + W, data=dat1)
y_stage1 <- model.matrix(~T + V + W, data=dat1)
lm_stage1 <- lm.fit(x=x_stage1, y=y_stage1)
dat1[, c("stage1_T", "stage1_V", "stage1_W", "stage1_1")] <-
  lm_stage1$fitted.values[, c("T0", "V1", "W1", "(Intercept)")]

# estimate stage 2
stage2 <- glm(
  evt ~ stage1_1 + stage1_T + stage1_V + stage1_W - 1,
  data = dat1,
  family = poisson(link="log")
)

res <- emmeans(stage2, ~ stage1_T, at=list(stage1_T=c(0,1))) |>
  pairs(type="response") |>
  summary(null=log(1-0.3), side="<", infer=TRUE, level=0.95)

# calculate VE as 1-risk ratio
1-res$ratio
1-res$asymp.UCL
res$p.value


# principal score weighting -----------------------------------------------

# prepare data
dat1 <- dat |>
  within({
    C <- factor(C, levels=c("0", "1"))
    trt <- factor(trt, levels=c("1", "0"))
  })

# propensity score model
mod_ps <- glm(C ~ V + W, subset = (trt==1), data=dat1, family=binomial())
# get odds from predicted values on the link (log-odds) scale
odds <- predict(mod_ps, newdata = dat1, type="link") |>
  exp()

# calculate weights
dat1 <- dat1 |>
  dplyr::mutate(
    weight = dplyr::case_when(
      ((trt==1) & (C=="1")) ~ 1,
      ((trt==0) & (C=="0")) ~ 0,
      ((trt==1) & (C=="0")) ~ 0,
      ((trt==0) & (C=="1")) ~ odds
    )
  )


# calculate outcomes model, weighted by score
# version 1, with V and W in the outcomes model
outcome_mod_1 <- glm(evt ~ trt + V + W, weights=weight, family=binomial(), data = dat1)

res <- emmeans(outcome_mod_1, ~ trt) |>
  regrid("log") |>
  contrast(method="pairwise", type="response") |>
  summary(null=log(1-0.3), side="<", infer=TRUE, level=0.95)

# calculate VE as 1-risk ratio
1-res$ratio
1-res$asymp.UCL
res$p.value

# calculate outcomes model, weighted by score
# version 2, without V and W in the outcomes model
outcome_mod_2 <- glm(evt ~ trt, weights=weight, family=binomial(), data = dat1)

res <- emmeans(outcome_mod_2, ~ trt) |>
  regrid("log") |>
  contrast(method="pairwise", type="response") |>
  summary(null=log(1-0.3), side="<", infer=TRUE, level=0.95)

# calculate VE as 1-risk ratio
1-res$ratio
1-res$asymp.UCL
res$p.value


# per protocol analysis ---------------------------------------------------

dat1 <- dat |>
  within({
    # recode as factor
    trt <- factor(trt, levels = c("1", "0"))
  })

mod_ve <- dat1 |>
  subset(C==1) |>
  glm(evt ~ trt + V + W, data = _, family=poisson(link="log"))

res <- emmeans(mod_ve, ~ trt) |>
  pairs(type="response") |>
  summary(null=log(1-0.3), side="<", infer=TRUE, level=0.95)

# calculate VE as 1-risk ratio
1-res$ratio
1-res$asymp.UCL
res$p.value
