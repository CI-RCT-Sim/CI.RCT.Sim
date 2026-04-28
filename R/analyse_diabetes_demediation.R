#' Create Analysis Function for De-mediation
#'
#' @param separate Default is TRUE which means that the imputation will be
#' performed separately for the treatment and control groups. If FALSE, the
#' imputation will be performed on the combined data set.
#'
#' @return an analyse function that can be used in runSimulation
#' @export
#'
#' @importFrom stats lm glm as.formula binomial coef predict var quantile
#' @importFrom mice mice make.method make.predictorMatrix complete rbind pool.scalar
#' @importFrom dplyr cur_column across matches mutate filter case_when select
#' @importFrom logistf logistf
#' @importFrom tidyselect starts_with
#'
#' @details
#' This function implements the de-mediation approach described in the paper
#' by Bartlett et al. (2024) for handling rescue medication in diabetes trials.
#' The function first imputes the missing rescue medication data using the
#' `mice` package, and then performs a de-mediation analysis by fitting a series
#' of models to estimate the effect of treatment on the outcome while accounting
#' for the potential mediation through rescue medication use.
#' The final estimate of the treatment effect is obtained by pooling the results
#' from the multiple imputations using Rubin's rules.
#'
#' @examples
#' \donttest{
#' setting <- diabetes_scenario()[1, ] |> diabetes_scenario_set_truevalues()
#' dat <- generate_diabetes(setting)
#' analyse_diabetes_demediation()(setting, dat)
#' }
analyse_diabetes_demediation <- function(separate = TRUE) {
  function(condition, dat, fixed_objects = NULL) {
    # Convert the logical of receiving rescue at any point in to a longitudinal measurement in wide format
    daat <- dat |> mutate(
      across(
        starts_with("R") & !R0 & !rescue_start,
        ~ ifelse(is.na(.) & rescue_start <= as.numeric(gsub("R", "", cur_column())), 1, .),
        .names = "R{gsub('R', '', .col)}"
      )
    )
    Rcols <- grep("^R[0-9]+$", names(daat), value = TRUE)

    if (separate) {
      dat1 <- daat |> filter(trt == 1)
      pred1 <- make.predictorMatrix(dat1)
      pred1[Rcols, ] <- 0
      pred1["m_start", ] <- 0
      pred1[, c("id", "trt", "m_start", Rcols)] <- 0
      meth1 <- make.method(dat1)
      meth1[Rcols] <- ""
      dat0 <- daat |> filter(trt == 0)
      pred0 <- make.predictorMatrix(dat0)
      pred0[Rcols, ] <- 0
      pred0["m_start", ] <- 0
      pred0[, c("id", "trt", "m_start", Rcols)] <- 0
      meth0 <- make.method(dat0)
      meth0[Rcols] <- ""
      dats <- rbind(
        mice(dat1, m = 5, printFlag = FALSE, predictorMatrix = pred1, method = meth1, ridge = 1e-5, remove.collinear = FALSE),
        mice(dat0, m = 5, printFlag = FALSE, predictorMatrix = pred0, method = meth0, ridge = 1e-5, remove.collinear = FALSE)
      )
    } else {
      pred <- make.predictorMatrix(daat)
      pred[Rcols, ] <- 0
      pred["m_start", ] <- 0
      pred[, c("id", "m_start", Rcols)] <- 0
      meth <- make.method(dat)
      meth[Rcols] <- ""
      dats <- mice(daat, m = 5, printFlag = FALSE, predictorMatrix = pred, method = meth, ridge = 1e-5, remove.collinear = FALSE)
    }

    analysis <- function(dataa) {
      dat_comp <- dataa |> mutate(
        across(
          matches("^y[0-9]+$") & !y0,
          ~ .x - y0,
          .names = "yc{gsub('y', '', .col)}"
        )
      )

      Rcols <- grep("^R[0-9]+$", names(dat_comp), value = TRUE)
      visits <- as.numeric(sub("R", "", Rcols))

      Rmat <- as.matrix(dat_comp[Rcols])
      rs <- dat_comp$rescue_start

      # NA & visit < rescue_start -> 0
      mask0 <- is.na(Rmat) & outer(rs, visits, `>`)

      # NA & visit >= rescue_start -> 1
      mask1 <- is.na(Rmat) & outer(rs, visits, `<=`)

      Rmat[mask0] <- 0
      Rmat[mask1] <- 1

      dat_comp[Rcols] <- Rmat

      dat_comp[, paste0("j", condition$k - 1)] <- dat_comp[, paste0("yc", condition$k)]
      for (k in 1:(condition$k - 1)) {
        if (sum(dat_comp[, paste0("R", condition$k - k)], na.rm = TRUE) == 0) {
          dat_comp[, paste0("j", condition$k - k - 1)] <- dat_comp[, paste0("j", condition$k - k)]
          next
        }

        # Fit a model to predict the probability of receiving rescue medication at visit 12 - k
        mod <- logistf(
          as.formula(paste0("R", condition$k - k, " ~ trt + age + y0", paste0("+ y", 1:(condition$k - k), collapse = " "))),
          data = dat_comp,
          pl = FALSE,
          control = logistf::logistf.control(maxit = 2000, maxstep = 0.5)
        )
        dat_comp[, paste0("pred_R", condition$k - k)] <- predict(mod, type = "response")

        # Subset the data
        daats <- dat_comp |>
          select(
            trt,
            age,
            y0,
            paste0("R", (condition$k - k):1),
            paste0("yc", (condition$k - k):1),
            paste0("pred_R", condition$k - k),
            paste0("j", condition$k - k)
          )

        # Fit a model for the outcome
        model <- lm(
          as.formula(
            paste0("j", condition$k - k, "~ .")
          ),
          data = daats
        )

        dat_comp[, paste0("j", condition$k - k - 1)] <-
          dat_comp[, paste0("j", condition$k - k)] -
          case_when(
            !is.na(coef(model)[paste0("R", condition$k - k)]) ~ coef(model)[paste0("R", condition$k - k)],
            TRUE ~ 0
          ) *
            dat_comp[, paste0("R", condition$k - k)]
      }

      final.model <- lm(j0 ~ trt + y0, data = dat_comp)
      c(
        p = summary(final.model)$coefficients["trt", "Pr(>|t|)"],
        coef = coef(final.model)["trt"],
        se = summary(final.model)$coefficients["trt", "Std. Error"]
      )
    }
    effect <- rep(NA, dats$m)
    effect.var <- rep(NA, dats$m)
    for (i in 1:dats$m) {
      dat <- complete(dats, i)
      res <- analysis(dat)
      effect[i] <- res["coef.trt"]
      effect.var[i] <- res["se"]^2
    }
    end_res <- pool.scalar(effect, effect.var)
    ci <- c(
      end_res$qbar - 1.96 * sqrt(end_res$t),
      end_res$qbar + 1.96 * sqrt(end_res$t)
    )

    list(
      ci_lower = ci[1],
      ci_upper = ci[2],
      coef = end_res$qbar,
      sd = sqrt(end_res$t)
    )
  }
}
