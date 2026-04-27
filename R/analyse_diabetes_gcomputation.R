#' Analyse data set with G-computation
#'
#' @return A function that, when called with `condition` and `dat`, returns a list with:
#' * `coef`      estimated difference in mean change in HbA1c between treatment groups
#' * `se`        standard error for coef
#' * `ci_lower`  lower bound of 95% confidence interval for coef
#' * `ci_upper`  upper bound of 95% confidence interval for coef
#'
#' @export
#'
#' @importFrom gfoRmula gformula_continuous_eof lagged static carry_forward
#' @importFrom dplyr mutate arrange group_by lead
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr `%>%`
#' @importFrom data.table as.data.table
#' @importFrom logistf logistf
#'
#' @details
#' This function implements G-computation using `gformula_continuous_eof` to estimate the
#' average treatment effect (ATE) on the **change in HbA1c** under a hypothetical scenario
#' where rescue medication is not available, regardless of treatment discontinuation.
#'
#' The estimand is:
#' - **Endpoint**: Change in HbA1c at final visit k
#' - **Treatments**: Experimental vs. control
#' - **Summary measure**: Difference in means
#' - **Population**: People with type 2 diabetes
#' - **Intercurrent events**: Use of rescue medication and treatment discontinuation
#' - **Intercurrent event strategy**:
#'   - *Hypothetical strategy*: Assume rescue medication was not available (i.e., rescue = 0 always)
#'   - *Treatment policy strategy*: Ignore treatment discontinuation (i.e., continue treatment in counterfactual)
#'
#' The joint distribution of age, HbA1c, and rescue medication is modeled over time:
#' - change HbA1c ~ lag1(HbA1c) + age + trt + R
#' - Rescue ~ change HbA1c + age
#'
#' A restriction is applied: once rescue is initiated, it remains on (to reflect continuous use).
#'
#' The intervention sets rescue = 0 for all patients and all time points (no rescue allowed).
#' Treatment discontinuation is ignored (treatment policy).
#'
#' Bootstrap (200 samples, percentile method) is used to compute 95% confidence intervals
#' for the mean difference in change in HbA1c between groups.
#'
#' The null hypothesis of no treatment effect is rejected if the 95% CI does not include 0.
#' No formal p-value is calculated.
#'
#' @examples
#' \donttest{
#' Design <- diabetes_scenario()[1, ] |>
#'   diabetes_scenario_set_truevalues()
#'
#' dat <- generate_diabetes(Design)
#'
#' analyse_diabetes_gcomputation()(Design, dat)
#' }
analyse_diabetes_gcomputation <- function() {
  function(Design, dat, fixed_objects = NULL) {
    k <- Design$k # number of last visit
    setup <- Design$setup # determines whether rescue medication is switched to (setup = 0) or put on top of active treatment (setup = 1)

    # reformat dat to long format with outcome column 'y' for the change in HbA1c at each visit
    dat_long <- pivot_longer(dat,
      cols = matches("^[yR]\\d+$"),
      names_to = c(".value", "visit"),
      names_pattern = "([yR])(\\d+)"
    ) %>%
      mutate(
        hba1c = y,
        visit = as.numeric(sub("y", "", visit)),
        R = ifelse(visit == 0, 0, R) # set rescue at baseline to 0
      ) %>%
      arrange(id, visit) %>%
      group_by(id) %>% # make sure table is grouped by id and ordered by visit
      mutate(
        hba1c_0 = hba1c[visit == 0], # HbA1c at baseline
        y = hba1c - hba1c_0, # HbA1c change
        # want to fit models on data up to time k-1, then simulate forward to predict the outcome at time k:
        # create a new column at second-to-last timepoint that holds y at last time point
        # i.e. final outcome
        # to preserve this value while deleting the last row for the model estimation
        y_k = ifelse(visit == k - 1, dplyr::lead(y), NA)
      ) # |>
    # dplyr::select(-R)

    # Remove final visit i.e. visit k
    dat_long <- dat_long[dat_long$visit != k, ]

    # Under assumption 1 (setup = 0)
    # when rescue is started
    # treatment group should switch to control, and control should stay in control
    if (setup == 0) {
      dat_long <- dat_long %>% mutate(trt = replace(trt, R == 1, 0))
    }

    # === Firth Correction: Replace glm for binomial models with logistf ===
    # This ensures robust fitting when quasi-separation occurs in rescue medication models
    firth_glm <- function(formula, family = gaussian(), data = NULL, ...) {
      if (inherits(family, "binomial")) {
        # Check for unsupported arguments
        unsupported <- c("offset", "weights", "subset", "na.action")
        if (any(sapply(unsupported, function(x) !is.null(get(x, envir = list(...)))))) {
          warning("logistf does not support offset, weights, subset, or na.action. Using default.")
        }
        # Fit with logistf (Firth's bias-reduced logistic regression)
        model <- tryCatch({
          logistf(formula, data = data, family = binomial, ...)
        }, error = function(e) {
          stop("Firth correction failed. Consider checking for separation or data issues.")
        })
        # Extract and return mimicking glm output
        list(
          coefficients = coef(model),
          se = sqrt(diag(vcov(model))),
          deviance = model$deviance,
          df.residual = model$df.residual,
          fitted.values = plogis(model$linear.predictors),
          model = model,
          call = match.call(),
          terms = terms(formula),
          na.action = attr(data, "na.action"),
          xlevels = model$xlevels,
          y = model$y,
          x = model$x,
          linear.predictors = model$linear.predictors,
          residuals = model$residuals,
          weights = model$weights,
          offset = model$offset,
          prior.weights = model$prior.weights,
          family = family,
          contrasts = model$contrasts,
          fitted.values = plogis(model$linear.predictors),
          ...
        )
      } else {
        # For non-binomial models, use standard glm
        return(glm(formula, family = family, data = data, ...))
      }
    }

    # === Run g-computation with Firth correction only for logistic models ===
    # Use local() to isolate the glm override
    g.model <- local({
      # Temporarily replace glm with safe_glm
      assign("glm", firth_glm, envir = .GlobalEnv)
      # Now call gformula_continuous_eof — it will use Firth for binomial models
      gformula_continuous_eof(
        obs_data = as.data.table(dat_long),
        id = "id",
        time_name = "visit",
        covnames = c("y", "trt", "R"),
        outcome_name = "y_k",
        covtypes = c("normal", "binary", "binary"),
        covparams = list(covmodels = c(
          y ~ trt + R + lag1_y + age,
          trt ~ 1,
          R ~ y + age
        )),
        ymodel = y_k ~ trt + R + lag1_y + age,
        intvars = list(c("trt", "R"), c("trt", "R")),
        interventions = list(
          list(c(static, rep(1, k)), c(static, rep(0, k))),
          list(c(static, rep(0, k)), c(static, rep(0, k)))
        ),
        int_descript = c("treatment no rescue", "control no rescue"),
        restrictions = list(c("R", "lag1_R != 1", carry_forward)),
        ref_int = 2,
        histvars = list("y", "R"),
        histories = list(lagged, lagged),
        basecovs = c("age", "hba1c_0"),
        nsamples = 200,
        show_progress = FALSE,
        #seed = 10,
        ci_method = "percentile"
      )
    })

    # mean of difference in mean change over bootstrap samples
    # summary(g.model)
    coef <- g.model$result$`Mean difference`[2] # mean difference between treatments (intervention - control) at visit k
    se <- g.model$result$`MD SE`[2] # se for mean difference
    ci_lower <- g.model[["result"]][["MD lower 95% CI"]][2] # 95% lower CI via percentile method
    ci_upper <- g.model[["result"]][["MD upper 95% CI"]][2] # 95% lower CI via percentile method

    list(
      coef = coef,
      se = se,
      ci_lower = ci_lower,
      ci_upper = ci_upper
    )
  }
}
