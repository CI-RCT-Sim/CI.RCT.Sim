#' Create Analyse Functions for ...
#'
#' @param X input can be used to pass parameters to the analyse function
#'
#' @return an analyse function that can be used in runSimulation
#' @export
#'
#' @importFrom stats lm glm as.formula binomial coef predict var quantile
#' @importFrom mice mice make.method make.predictorMatrix complete rbind pool.scalar
#' @importFrom logistf logistf
#'
#' @examples
#' setting <- assumptions_diabetes_rescue()[1, ] |> true_summary_statistics_diabetes_rescue()
#' dat <- generate_diabetes_rescue(setting)
#' analyse_diabetes_demediation()(setting, dat)
#'
analyse_diabetes_demediation <- function(X) {
  function(condition, dat, fixed_objects = NULL) {
    # Convert the logical of receiving rescue at any point in to a longitudinal measurement in wide format
    dat <- dat |> dplyr::mutate(rescue_start = ifelse(is.na(rescue_start), condition$k + 2, rescue_start),R0 = 0)

    dat1 <- dat |> dplyr::filter(trt == 1)
    pred1 <- mice::make.predictorMatrix(dat1)
    pred1[, c("id", "rescue_start", "trt")] <- 0
    dat0 <- dat |> dplyr::filter(trt == 0)
    pred0 <- mice::make.predictorMatrix(dat0)
    pred0[, c("id", "rescue_start", "trt")] <- 0
    dats <- mice::rbind(
      mice::mice(dat1, m = 5, print = FALSE, predictorMatrix = pred1),
      mice::mice(dat0, m = 5, print = FALSE, predictorMatrix = pred0)
    )

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
        if (sum(dat_comp[, paste0("R", 12 - k)], na.rm = TRUE) == 0) {
          dat_comp[, paste0("j", 12 - k - 1)] <- dat_comp[, paste0("j", 12 - k)]
          next
        }
        # Fit a model to predict the probability of receiving rescue medication at visit 12 - k
        # browser()
        # mod <- glm(
        #   as.formula(paste0("R", 12 - k, "~ trt + age + y0", paste0("+ yc", 1:(12 - k), collapse = " "))),
        #   data = dat_comp,
        #   family = binomial(link = "probit"),
        #   method = "brglmFit"
        # )
        # mod <- logistf(
        #   as.formula(paste0("R", 12 - k, "~ trt + age + y0", paste0("+ yc", 1:(12 - k), collapse = " "))),
        #   data = dat_comp
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
      # browser()
      c(
        p = summary(final.model)$coefficients["trt", "Pr(>|t|)"],
        coef = coef(final.model)["trt"],
        se = summary(final.model)$coefficients["trt", "Std. Error"]
      )
    }
    effect <- rep(NA, dats$m)
    effect.var <- rep(NA, dats$m)
    cil <- rep(NA, dats$m)
    ciu <- rep(NA, dats$m)
    # res <- vector("list", dats$m)
    # browser()
    for (i in 1:dats$m) {
      dat <- mice::complete(dats, i)
      res <- boot::boot(dat, analysis, R = 500)
      effect[i] <- res$t0[2]
      effect.var[i] <- var(res$t[, 2])
      cil[i] <- quantile(res$t[, 2], probs = c(0.025))
      ciu[i] <- quantile(res$t[, 2], probs = c(0.975))
    }
    end_res <- mice::pool.scalar(effect, effect.var)

    list(
      effect,
      effect.var,
      cil,
      ciu,
      coefs = end_res$qhat,
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
