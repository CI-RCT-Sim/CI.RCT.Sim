#' Create Analyse Functions for ...
#'
#' @param X input can be used to pass parameters to the analyse function
#'
#' @return an analyse function that can be used in runSimulation
#' @export
#'
#' @importFrom stats lm glm as.formula binomial coef predict var quantile
#' @importFrom mice mice make.method make.predictorMatrix complete rbind pool.scalar
#' @importFrom dplyr cur_column across matches mutate filter case_when select
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
    daat <- dat |> mutate(
      rescue_start = ifelse(is.na(rescue_start), condition$k + 2, rescue_start),
      R0 = 0, .before = R1,
      across(
        starts_with("R") & !R0 & !rescue_start,
        ~ ifelse(is.na(.) & rescue_start <= as.numeric(gsub("R", "", cur_column())), 1, .),
        .names = "R{gsub('R', '', .col)}"
      )
    )

    dat1 <- daat |> filter(trt == 1)
    pred1 <- make.predictorMatrix(dat1)
    pred1[, c("id", "rescue_start", "trt")] <- 0
    meth1 <- make.method(dat1)
    dat0 <- daat |> filter(trt == 0)
    pred0 <- make.predictorMatrix(dat0)
    pred0[, c("id", "rescue_start", "trt")] <- 0
    meth0 <- make.method(dat0)
    dats <- rbind(
      mice(dat1, m = 10, print = FALSE, predictorMatrix = pred1, method = meth1, ridge = 1e-5, remove.collinear = FALSE),
      mice(dat0, m = 10, print = FALSE, predictorMatrix = pred0, method = meth0, ridge = 1e-5, remove.collinear = FALSE)
    )
    # browser()
    analysis <- function(dataa) {
      dat_comp <- dataa |> mutate(
        across(
          matches("^y[0-9]+$") & !y0,
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
        mod <- logistf(
          as.formula(paste0("R", 12 - k, "~ trt + age + y0", paste0("+ yc", 1:(12 - k), collapse = " "))),
          data = dat_comp,
          pl = FALSE
        )
        browser(expr = length(predict(mod, type = "response")) != nrow(dat_comp))
        dat_comp[, paste0("pred_R", 12 - k)] <- predict(mod, type = "response")

        # Subset the data
        daat <- dat_comp |>
          select(
            trt,
            age,
            y0,
            paste0("R", (12 - k):1),
            paste0("yc", (12 - k):1),
            paste0("pred_R", 12 - k),
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
          case_when(
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
    ests <- list()
    # res <- vector("list", dats$m)
    # browser()
    for (i in 1:dats$m) {
      dat <- complete(dats, i)
      dat[is.na(dat)] <- 0
      # browser(expr = sum(is.na(dat)) != 0)
      res <- boot::boot(dat, analysis, R = 500)
      res <- analysis(dat)
      # effect[i] <- res$t0[2]
      # effect.var[i] <- var(res$t[, 2])
      # ests[[i]] <- res$t[, 2]
      effect[i] <- res["coef.trt"]
      effect.var[i] <- res["se"]^2
    }
    end_res <- pool.scalar(effect, effect.var)
    ci <- c(
      end_res$qbar - 1.96 * sqrt(end_res$t),
      end_res$qbar + 1.96 * sqrt(end_res$t)
    )
    # browser()
    list(
      # effect,
      # effect.var,
      ci,
      # cil = quantile(unlist(ests), 0.025),
      # ciu = quantile(unlist(ests), 0.975),
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
