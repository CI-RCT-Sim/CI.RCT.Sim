#' Create Analyse Functions for ...
#'
#' @param X
#'
#' @return an analyse function that can be used in runSimulation
#' @export
#'
#' @importFrom stats lm glm as.formula
#'
#' @examples
#' analyse_demediation()(assumptions_diabetes_rescue()[1, ], generate_diabetes_rescue(assumptions_diabetes_rescue()[1, ]))
#'
analyse_demediation <- function(X) {
  function(condition, dat, fixed_objects = NULL) {
    # Convert the logical of receiving rescue at any point in to a longitudinal measurement in wide format
    dat <- dat |> dplyr::mutate(rescue_start = ifelse(is.na(rescue_start), condition$k[1] + 2, rescue_start))
    pred <- mice::make.predictorMatrix(dat)
    pred[upper.tri(pred)] <- 0
    dats <- mice::mice(dat,
      m = 2, print = FALSE, seed = 2026,
      formulas = mice::make.formulas(dat, blocks = mice::make.blocks(dat), predictorMatrix = pred)
    )
    # t <- lapply(1:dats$m, function(i) {
    #   dat <- mice::complete(dats, i)
    #   dat_comp <- dat |> dplyr::mutate(
    #     across(
    #       matches("^y[0-9]+$") & !y0,
    #       ~ .x - y0,
    #       .names = "yc{gsub('y', '', .col)}"
    #     )
    #   )
    #   for (j in 1:condition$k[1]) {
    #     dat_comp[[paste0("R", j)]] <- ifelse(dat_comp$rescue_start <= j, 1L, 0L)
    #   }
    #
    #   # Step 1: Set Y at T_10 as R_t
    #   M.t <- dat_comp$y12
    #
    #   k <- 1
    #   for (k in 1:11) {
    #     data <- dat_comp |>
    #       # Only treatment, baseline characteristics
    #       # And repeated HbA1c measurements
    #       # ICE indicator as treatment
    #       dplyr::select(
    #         trt,
    #         "y0",
    #         paste0("R", 1:(12 - k)),
    #         paste0("y", 1:(12 - k))
    #       ) |>
    #       as.data.frame()
    #
    #     # note that we here use ICE indicators which are defined as E's in the paper,
    #     # whereas our description of G-estimation there refers to M's. However,
    #     # due to the relation between the M's and E's, it makes no difference to
    #     # the estimator whether E's or M's are used
    #
    #     g.est.model <-
    #       lm(M.t ~ .,
    #         data = data
    #       )
    #
    #     mediator <-
    #       as.data.frame(data[, paste0("R", 12 - k)])
    #     # in some bootstrap samples sometimes the coefficient required is not estimable
    #     # in these cases we set the coefficient to zero, and so R.t is no longer removing
    #     # the effect of the corresponding ICE variable
    #     psi <- ifelse(!is.na(coef(g.est.model)[paste0("R", 12 - k)]),
    #       coef(g.est.model)[paste0("R", 12 - k)], 0
    #     )
    #
    #     M.t <- as.vector(M.t -
    #       psi * mediator)[[1]]
    #   }
    #   g.model <-
    #     lm(M.t ~ data$trt + data$y0)
    #   g.model
    #   list(
    #     p = summary(g.model)$coefficients["data$trt", "Pr(>|t|)"],
    #     coef = coef(g.model)["data$trt"],
    #     ci_lower = confint(g.model)["data$trt", 1],
    #     ci_upper = confint(g.model)["data$trt", 2],
    #     iteration = i
    #   )
    # })
    # browser()

    t <- lapply(1:dats$m, function(i) {
      dat <- mice::complete(dats, i)
      dat_comp <- dat |> dplyr::mutate(
        across(
          matches("^y[0-9]+$") & !y0,
          ~ .x - y0,
          .names = "yc{gsub('y', '', .col)}"
        )
      )
      for (j in 1:condition$k[1]) {
        dat_comp[[paste0("R", j)]] <- ifelse(dat_comp$rescue_start <= j, 1L, 0L)
      }

      dat_comp[, "j11"] <- dat_comp$yc12
      for (k in 1:11) {
        mod <- glm(as.formula(paste0("R", 12 - k, "~ trt + y0")), data = dat_comp, family = binomial(link = "probit"))
        dat_comp[, paste0("pred_R", 12 - k)] <- predict(mod, type = "response")
        model <- lm(as.formula(paste0("j", 12 - k, "~ trt + y0 + pred_R", 12 - k, " + R", 12 - k)), data = dat_comp)
        dat_comp[, paste0("j", 12 - k - 1)] <- dat_comp[, paste0("j", 12 - k)] - coef(model)[paste0("R", 12 - k)] * dat_comp[, paste0("R", 12 - k)]
      }

      final.model <- lm(j0 ~ trt + y0, data = dat_comp)
      list(
        p = summary(final.model)$coefficients["trt", "Pr(>|t|)"],
        coef = coef(final.model)["trt"],
        ci_lower = confint(final.model)["trt", 1],
        ci_upper = confint(final.model)["trt", 2]
        # iteration = i
      )
    })

    # analysis <- function(dat) {
    #   # dat <- mice::complete(dats, i)
    #   dat_comp <- dat |> dplyr::mutate(
    #     across(
    #       matches("^y[0-9]+$") & !y0,
    #       ~ .x - y0,
    #       .names = "yc{gsub('y', '', .col)}"
    #     )
    #   )
    #   for (j in 1:12) {
    #     dat_comp[[paste0("R", j)]] <- ifelse(dat_comp$rescue_start <= j, 1L, 0L)
    #   }
    #
    #   dat_comp[, "j11"] <- dat_comp$yc12
    #   for (k in 1:11) {
    #     mod <- glm(as.formula(paste0("R", 12 - k, "~ trt + y0")), data = dat_comp, family = binomial(link = "probit"))
    #     dat_comp[, paste0("pred_R", 12 - k)] <- predict(mod, type = "response")
    #     model <- lm(as.formula(paste0("j", 12 - k, "~ trt + y0 + pred_R", 12 - k, " + R", 12 - k)), data = dat_comp)
    #     dat_comp[, paste0("j", 12 - k - 1)] <- dat_comp[, paste0("j", 12 - k)] - coef(model)[paste0("R", 12 - k)] * dat_comp[, paste0("R", 12 - k)]
    #   }
    #
    #   final.model <- lm(j0 ~ trt + y0, data = dat_comp)
    #   list(
    #     p = summary(final.model)$coefficients["trt", "Pr(>|t|)"],
    #     coef = coef(final.model)["trt"],
    #     ci_lower = confint(final.model)["trt", 1],
    #     ci_upper = confint(final.model)["trt", 2]
    #   )
    # }
    # t <- with(dats,analysis)

    # browser()
    as.list(dplyr::bind_rows(t)[1, ])
  }
}


#' Summarise Output from Analyse Functions for ...
#'
#' @param X
#'
#' @describeIn analyse_demediation Summarise Output from Analyse X
#'
#' @return
#' Returns a function with the arguments:
#'  * condition
#'  * results
#'  * fixed objects
#'
#' that can be passed to create_summarise_function or to
#' SimDesign::runSimulation and that returns a `data.frame` with the columns
#'  * `Y` ...
#'  * ...
#'
#' @export
#'
#' @examples
#' summarise_demediation("tell")
summarise_demediation <- function(name = NULL) {
  # res <- function(condition, results, fixed_objects = NULL) {
  res <- data.frame(
    "Y" = NA_real_
  )
  # }

  attr(res, "name") <- name

  res
}
