#' Analyse data set with G-computation
#'
#' @return A function that, when called with `condition` and `dat`, returns a list with:
#' * `coef`      estimated difference in mean change in HbA1c between treatment groups
#' * `se`        standard error for coef
#' * `p`         p-value for coef
#' * `ci_lower`  lower bound of 95% confidence interval for coef
#' * `ci_upper`  upper bound of 95% confidence interval for coef
#'
#' @export
#'
#' @importFrom dplyr mutate arrange group_by ungroup lead lag filter select
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr `%>%`
#'
#' @details
#' This function implements G-computation to estimate the
#' average treatment effect (ATE) on the **change in HbA1c** under a hypothetical scenario
#' where rescue medication is not available, regardless of treatment discontinuation.
#'
#' Two scenarios are considered depending on `setup`:
#' - `setup = 1`: the effect of rescue goes on top of the effect of the assigned treatment
#' - `setup = 0`: once rescue starts, treatment is deterministically recoded to control
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
# Data are reshaped to long format. Baseline HbA1c, change from baseline,
# lagged HbA1c, and the final outcome y_k are constructed, and the final visit
# row is removed so that visit k - 1 contains the full covariate history needed
# to model the end-of-follow-up outcome.
#
# The fitted models are:
#   hba1c ~ trt + R + lag1_hba1c + age
#   y_k   ~ trt + R + hba1c + age + hba1c_0
#
# Counterfactual trajectories are then simulated separately for each treatment arm,
# starting from the observed baseline covariate distribution and imposing R = 0
# at all post-baseline visits. HbA1c is simulated forward under the intervention,
# and the final change in HbA1c is simulated from the outcome model.
# The treatment effect is the difference in mean simulated final change between arms.
#
# A separate rescue model is not used in the simulation step because the intervention
# is deterministic: rescue is fixed to 0 at all times. Thus, for this specific
# no-rescue estimand, the post-intervention distribution can be simulated directly
# from the HbA1c and final outcome models evaluated at R = 0. A rescue model would
# be needed for a more general g-formula implementation of the natural course or
# alternative rescue strategies.
#
# Uncertainty is quantified by bootstrapping (200 samples, percentile method) the entire procedure:
# resample subjects, rebuild the long data, refit the models, resimulate the
# no-rescue counterfactual outcomes, and recompute the treatment contrast.
#'
#' The null hypothesis of no treatment effect is rejected if the 95% CI does not include 0.
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
  function(condition, dat, fixed_objects = NULL) {
    k <- condition$k # number of last visit
    setup <- condition$setup # determines whether rescue medication is switched to (setup = 0) or put on top of active treatment (setup = 1)

    estimate_gformula <- function(dat) { # this function is run in each bootstrap iteration

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
          y_k = ifelse(visit == k - 1, dplyr::lead(y), NA),
          lag1_hba1c = dplyr::lag(hba1c),
          lag1_R = dplyr::lag(R)
        ) %>%
        ungroup()

      # Remove final visit i.e. visit k
      dat_long <- dat_long[dat_long$visit != k, ]

      # Under assumption 1 (setup = 0)
      # when rescue is started
      # treatment group should switch to control, and control should stay in control
      if (setup == 0) {
        dat_long <- dat_long %>% mutate(trt = replace(trt, R == 1, 0))
      }
      # Run g-computation with bootstrap using custom function
      # We simulate HbA1c values under the intervention (no rescue) and then
      # estimate the mean change in HbA1c at the final visit


      dat_followup <- dat_long %>%
        filter(visit > 0)

      # fit visit-specific hba1c models
      fit_hba1c <- vector("list", length = k - 1)

      for (v in 1:(k - 1)) {
        dat_visit <- dat_followup %>%
          filter(visit == v)

        fit_hba1c[[v]] <- stats::lm(
          hba1c ~ trt + R + lag1_hba1c + age,
          data = dat_visit
        )
      }

      # retain only second to last hba1c value
      dat_outcome <- dat_long %>%
        filter(visit == k - 1)

      # fit outcome model (change in hba1c)
      fit_ymodel <- stats::lm(
        y_k ~ trt + R + hba1c + age + hba1c_0,
        data = dat_outcome
      )

      baseline_dat <- dat_long %>%
        filter(visit == 0) %>%
        select(id, age, hba1c_0)

      simulate_mean <- function(trt_value) {
        sim_dat <- baseline_dat %>%
          mutate(
            trt = trt_value,
            hba1c = hba1c_0,
            R = 0
          )
        # iterate over each timepoint after baseline until k-1
        # put in values under intervention of no rescue instead of original observed data
        if (k > 1) {
          for (v in 1:(k - 1)) {
            newdat_hba1c <- data.frame(
              trt = sim_dat$trt,
              R = 0,
              lag1_hba1c = sim_dat$hba1c,
              age = sim_dat$age
            )
            mu_hba1c <- stats::predict(fit_hba1c[[v]], newdata = newdat_hba1c)
            sim_dat$hba1c <- stats::rnorm(nrow(sim_dat), mean = mu_hba1c, sd = summary(fit_hba1c[[v]])$sigma)
            sim_dat$R <- 0
          }
        }
        # based on predicted hypothetical follow-up values, simulate hypothetical outcome
        newdat_y <- data.frame(
          trt = sim_dat$trt,
          R = 0,
          hba1c = sim_dat$hba1c,
          age = sim_dat$age,
          hba1c_0 = sim_dat$hba1c_0
        )
        mu_y <- stats::predict(fit_ymodel, newdata = newdat_y)
        y_k_sim <- stats::rnorm(nrow(sim_dat), mean = mu_y, sd = summary(fit_ymodel)$sigma)

        mean(y_k_sim)
      }

      mean_trt <- simulate_mean(trt_value = 1)
      mean_ctrl <- simulate_mean(trt_value = 0)

      mean_trt - mean_ctrl
    }

    coef <- estimate_gformula(dat)

    # variance estimation using full bootstrap
    nsamples <- 200
    ids <- unique(dat$id)
    n_ids <- length(ids)
    boot_est <- numeric(nsamples)

    for (b in seq_len(nsamples)) {
      boot_ids <- sample(ids, size = n_ids, replace = TRUE)

      dat_boot <- do.call(rbind, lapply(seq_along(boot_ids), function(i) {
        dat_i <- dat[dat$id == boot_ids[i], ]
        dat_i$id <- i
        dat_i
      }))

      boot_est[b] <- estimate_gformula(dat_boot)
    }

    # mean of difference in mean change over bootstrap samples
    se <- stats::sd(boot_est) # se for mean difference
    ci_lower <- as.numeric(stats::quantile(boot_est, probs = 0.025)) # 95% lower CI via percentile method
    ci_upper <- as.numeric(stats::quantile(boot_est, probs = 0.975)) # 95% upper CI via percentile method
    p <- 2 * (1 - pnorm(abs(coef / se)))

    list(
      coef = coef,
      se = se,
      p = p,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
    )
  }
}
