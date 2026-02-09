#' Analyse Dataset with G-estimation (G-computation via gformula_continuous_eof)
#'
#' @param level confidence level for CI computation (default: 0.95)
#' @param alternative alternative hypothesis for the test (default: "two.sided")
#'   - "two.sided": tests H0: treatment effect = 0 vs H1: treatment effect ≠ 0
#'   - "one.sided": tests H0: treatment effect ≤ 0 vs H1: treatment effect > 0
#'   (Note: No formal p-value is computed; inference is based on CI)
#'
#' @return A function that, when called with `condition` and `dat`, returns a list with:
#' * `p`         NA (no formal p-value computed)
#' * `alternative` the alternative hypothesis used
#' * `coef`      estimated difference in mean change in HbA1c between treatment groups
#' * `ci_lower`  lower CI for mean change in HbA1c
#' * `ci_upper`  upper CI for mean change in HbA1c
#' * `CI_level`  confidence level used
#' * `N_pat`     number of patients in the dataset
#'
#' @export
#'
#' @importFrom gfoRmula gformula_continuous_eof
#' @importFrom dplyr filter lag arrange group_by
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr `%>%`
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
#' - HbA1c ~ lag1(HbA1c) + age + trt
#' - Rescue ~ HbA1c + age + trt
#'
#' A restriction is applied: once rescue is initiated, it remains on (to reflect continuous use).
#' The intervention sets rescue = 0 for all patients and all time points (no rescue allowed).
#' Treatment discontinuation is ignored (treatment policy).
#'
#' Bootstrap (500 samples, percentile method) is used to compute 95% confidence intervals
#' for the mean difference in change in HbA1c between groups.
#'
#' The null hypothesis of no treatment effect is rejected if the 95% CI does not include 0.
#' No formal p-value is calculated.
#'
#' @examples
#' Design <- assumptions_diabetes_rescue() |>
#'   true_summary_statistics_diabetes_rescue()
#'
#' condition <- Design[1, ]
#'
#' dat <- generate_diabetes_rescue(condition)
#'
#' analyse_gestimation <- analyse_gestimation()
#' analyse_gestimation(condition, dat)


analyse_gestimation <- function(level = 0.95, alternative = "two.sided") {
  stopifnot(alternative %in% c("two.sided", "one.sided"))
  stopifnot(level > 0 & level < 1)

  # What still remains to be done is
  # delete
  # add restriction to gestimation function
  # fix bootstrap (draw ids with replacement --> ids need to get new names otherwise the visits are messed up)
  # calculate mean over bootstrap samples
  # confidence interval
  # documentation
  # determine if the "alternative" statement can be incorporated

  function(condition, dat, fixed_objects = NULL) {
    k<-condition$k[1] #number of last visit

    #reformate dat to long format with one outcome column 'y' and a time variable 'visit'
    dat_long <- pivot_longer(dat, y0:paste0("y",k), names_to = "visit", values_to = "y")

    dat_long<- dat_long %>%
      mutate(visit = as.numeric(sub('hba1c', '', visit))) %>%
      mutate(rescue = ifelse(!is.na(rescue_start)&rescue_start<=visit, 1, 0)) %>% #new variable for rescue at visit j
      mutate(hba1c_0 = hba1c[visit == 0]) %>% # HbA1 at baseline
      mutate(y = hba1c - hba1c_0) %>% # HbA1c change
      arrange(id, visit) %>% group_by(id) %>% # make sure table is grouped by id and ordered by visit

    # Run g-computation with bootstrap
    # We simulate Hba1c values under the intervention (no rescue) and then
    # estimate the mean change in HbA1c at the final visit
#library('Hmisc')

    #create lists for bootstrap results
    coef_boot <- list()

    for (b in 1:500) { #bootstrap
      #draw bootstrap sample from patients
      obs <- unique(dat_long$id)
      obs_inds <- data.frame(id=sample(obs, length(obs), replace=TRUE))
      boot_data <- merge(obs_inds, dat_long, all.x=TRUE)

      gcomp_result<-list() #for each treatment
      for (i in c(0,1)){ # for each treatment group
        #gcomputation
        gcomp_result[i+1] <- gfoRmula::gformula_continuous_eof(
          obs_data = boot_data[boot_data$trt == i,],
          id = "id",
          time_name = "visit",
          covnames = c("y", "rescue"),
          covtypes = c("normal", "binary"),
          outcome_name = "y",
          basecovs = c("age"),
          histories = c(lagged),
          histvars = list(c('y')),
          covparams = list(covmodels = c(y ~ lag1_y + age, #for y
                                         rescue ~ y + age)), #for rescue
          #restrictions=list(c("rescue",condition, function, value used by function))
          ymodel = y ~ lag1_y + age,
          intervention1.rescue = list(static, rep(0, max(dat_long$visit) + 1)), # hypothetical: no rescue at any visit
          int_descript = c("hypothetical: no rescue was taken"),
          nsimul = length(unique(boot_data[boot_data$trt == i,]$id)),
          seed = 1234,
          nsamples = 0, #bootstrap will be conducted manually
          sim_data_b=TRUE,
          show_progress = TRUE
        )
        # Extract hypothetical final HbA1c and baseline HbA1c for each individual
        sim_data[i+1] <- gcomp_result[i+1]$sim_data$`hypothetical: no rescue was taken` %>% # intervention: no rescue
                          group_by(id) %>%
                          mutate(y_k = y[visit == k])%>%
                          mutate(y_0 = y[visit == 0])
        # Compute individual change
        sim_data[i+1] <-  sim_data[i+1] %>%
          mutate(y_change = y_k - y_0)

        # Compute mean change for each treatment group
        mean_change[i+1] <- mean(sim_data[i+1]$y_change, na.rm = TRUE)
      }

    # Difference in mean change between groups
    coef_boot[b] <- mean_change[1] - mean_change[0]
    }

    # mean of difference in mean change over bootstrap samples
    # confidence interval via

    list(
      p = NA,
      alternative = alternative,
      coef = mean_coef_boot,
      hr = NA,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      CI_level = level,
      N_pat = nrow(dat),
      N_evt = NA
    )
  }
}


