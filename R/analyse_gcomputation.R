#' Analyse Dataset with G-estimation (G-computation via gformula_continuous_eof)
#'
#' @return A function that, when called with `condition` and `dat`, returns a list with:
#' * `coef`      estimated difference in mean change in HbA1c between treatment groups
#' * `sd`       standard error for coef
#'
#' @export
#'
#' @importFrom gfoRmula gformula_continuous_eof
#' @importFrom dplyr arrange group_by
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr `%>%`
#' @importFrom data.table as.data.table
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
#' - change HbA1c ~ lag1(HbA1c) + age + trt
#' - Rescue ~ change HbA1c + age + trt
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


analyse_gestimation <- function() {

  # What still remains to be done is
  # add restriction to gestimation function
  # find out what the output is we are interested in
  # documentation

  function(condition, dat, fixed_objects = NULL) {
    k<-condition$k[1] #number of last visit

    #reformate dat to long format with outcome column 'y' for the change in HbA1c at each visit
    dat_long <- pivot_longer(dat, y0:paste0("y",k), names_to = "visit", values_to = "hba1c")

    dat_long<- dat_long %>%
      mutate(visit = as.numeric(sub('y', '', visit))) %>%
      mutate(rescue = ifelse(!is.na(rescue_start)&rescue_start<=visit, 1, 0))%>% #new variable for rescue at visit j
      arrange(id, visit) %>% group_by(id) %>% # make sure table is grouped by id and ordered by visit
      mutate(hba1c_0 = hba1c[visit == 0]) %>% # HbA1 at baseline
      mutate(y = hba1c - hba1c_0) # HbA1c change

    # want to fit models on data up to time k-1, then simulate forward to predict the outcome at time k:
    # create a new column at second-to-last timepoint that holds y at last time point
    # i.e. final outcome
    # to preserve this value while deleting the last row for the model estimation
    dat_long <- dat_long %>%
      group_by(id) %>%
      mutate(y_k = ifelse(visit == k - 1, lag(y), NA)) %>%
      ungroup()

    # Remove final visit i.e. visit k
    dat_long <- dat_long[dat_long$visit!=k,]


    # Run g-computation with bootstrap
    # We simulate Hba1c values under the intervention (no rescue) and then
    # estimate the mean change in HbA1c at the final visit

        # Parameters for g-formula function
        id <- "id"
        obs_data <- as.data.table(dat_long)
        time_name <- "visit"
        time_points <- k + 1 # number of time-points (+1 because baseline is included)
        covnames <- c("y", "trt", "rescue")
        outcome_name <- "y_k"
        covtypes <- c("normal","binary", "binary")
        histories <- c(lagged)
        histvars <- list("y", "rescue")
        basecovs <- c("age")
        covparams <- list(covmodels = c(y ~  trt + rescue + lag1_y + age, # + rescue??? (not in the protocol but would make sense to me becasue recue is affected by y-1 and affects y in future)
                                        trt ~ 1,
                                        rescue ~ trt + y + age)) # correct to put treatment in here?
        ymodel <- y_k ~ trt + rescue + lag1_y + age
        intvars <- list(c('trt', 'rescue'),
                        c('trt', 'rescue'))
        interventions <- list(list(c(static, rep(1, time_points)), #treatment
                                   c(static, rep(0, time_points))), #no rescue
                              list(c(static, rep(0, time_points)), #no treatment
                                   c(static, rep(0, time_points)))) #no rescue
        int_descript <- c('Hypothetical treatment no rescue',
                          'Hypothetical control no rescue')
        restrictions <- list(c("rescue","lag1_rescue == 1","replace",1))
        nsamples <- 5 # 500

        g.model <- gformula_continuous_eof(
                                  obs_data = obs_data,
                                  id = id,
                                  time_name = time_name,
                                  covnames = covnames,
                                  outcome_name = outcome_name,
                                  covtypes = covtypes,
                                  covparams = covparams,
                                  ymodel = ymodel,
                                  intvars = intvars,
                                  interventions = interventions,
                                  int_descript = int_descript,
                                  restrictions=restrictions,
                                  ref_int = 1,
                                  histvars = histvars,
                                  histories = histories,
                                  basecovs = basecovs,
                                  nsamples = nsamples,
                                  show_progress = F,
                                  seed = 1)

    # mean of difference in mean change over bootstrap samples
    # standard error

    list(
      coef = NA,
      se = NA
    )
  }
}


