#' Create Analyse Functions for ...
#'
#' @param X input can be used to pass parameters to the analyse function
#'
#' @return an analyse function that can be used in runSimulation
#' @export
#'
#' @importFrom stats lm glm as.formula binomial coef predict var
#'
#' @examples
#' setting <- assumptions_diabetes_rescue()[1, ]
#' dat <- generate_diabetes_rescue(setting)
#' analyse_diabetes_demediation()(setting, dat)
#'
analyse_diabetes_demediation <- function(X) {
  function(condition, dat, fixed_objects = NULL) {
    # Convert the logical of receiving rescue at any point in to a longitudinal measurement in wide format
    dat <- dat |> dplyr::mutate(rescue_start = ifelse(is.na(rescue_start), condition$k + 2, rescue_start))
    dats <- mice::mice(dat, m = 5, print = FALSE)

    analysis <- function(dat, indicator) {
      dat_comp <- dat[indicator, ] |> dplyr::mutate(
        dplyr::across(
          dplyr::matches("^y[0-9]+$") & !y0,
          ~ .x - y0,
          .names = "yc{gsub('y', '', .col)}"
        )
      )

      dat_comp[, "j11"] <- dat_comp$yc12
      for (k in 1:11) {
        # Fit a model to predict the probability of receiving rescue medication at visit 12 - k
        # mod <- glm(
        #   as.formula(paste0("R", 12 - k, "~ trt + age + y0", paste0("+ yc", 1:(12 - k), collapse = " "))),
        #   data = dat_comp,
        #   family = binomial(link = "probit")
        # )
        # dat_comp[, paste0("pred_R", 12 - k)] <- predict(mod, type = "response")

        # Subset the data
        daat <- dat_comp |>
          dplyr::select(
            trt,
            age,
            y0,
            paste0("R", (12 - k):1),
            paste0("yc", (12 - k):1),
            # paste0("pred_R", 12 - k),
            paste0("j", 12 - k)
          )

        # Fit a model for the outcome
        model <- lm(
          as.formula(
            paste0("j", 12 - k, "~ .")
            # paste0("j", 12 - k, "~ trt + age + y0 + pred_R", 12 - k, " + R", 12 - k)
          ),
          data = daat
        )
        dat_comp[, paste0("j", 12 - k - 1)] <-
          dat_comp[, paste0("j", 12 - k)] -
          dplyr::case_when(
            !is.na(coef(model)[paste0("R", 12 - k)]) ~ coef(model)[paste0("R", 12 - k)],
            TRUE ~ 0
          ) *
            dat_comp[, paste0("R", 12 - k)]
      }


      final.model <- lm(j0 ~ trt + y0, data = dat_comp)
      # browser(expr = abs(final.model$coefficients["trt"]) > 1)
      c(
        p = summary(final.model)$coefficients["trt", "Pr(>|t|)"],
        coef = coef(final.model)["trt"]
      )
    }
    effect <- rep(NA, dats$m)
    effect.var <- rep(NA, dats$m)
    # browser()
    for (i in 1:dats$m) {
      dat <- mice::complete(dats, i)
      res <- boot::boot(dat, analysis, R = 500)
      effect[i] <- res$t0[2]
      effect.var[i] <- var(res$t[, 2])
    }
    end_res <- mice::pool.scalar(effect, effect.var)

    list(
      # effect,
      # effect.var,
      # coefs = end_res$qhat,
      coef = end_res$qbar,
      sd = sqrt(end_res$t)
    )
  }
}


#' Summarise Output from Analyse Functions for ...
#'
#' @param name also input used to name the summarise function
#'
#' @describeIn summarise_diabetes_demediation Summarise Output from Analyse X
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
#' summarise_diabetes_demediation("tell")
summarise_diabetes_demediation <- function(name = NULL) {
  # res <- function(condition, results, fixed_objects = NULL) {
  res <- data.frame(
    "Y" = NA_real_
  )
  # }

  attr(res, "name") <- name

  res
}
