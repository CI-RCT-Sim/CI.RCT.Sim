analyse_diabetes_rescue_mi <- function(m = 10,
                                       maxit = 20,
                                       ci_level = 0.95,
                                       seed = 123) {

  function(condition, dat, fixed_objects = NULL) {

    library(mice)
    library(dplyr)

    k <- condition$k

    vars_y <- paste0("y", 0:k)
    vars_R <- paste0("R", 1:(k - 1))

    vars_imp <- c(vars_y, vars_R, "age", "trt")

    # Methods
    meth <- make.method(dat[vars_imp])

    # HbA1c imputation
    meth[vars_y] <- "pmm"

    # Rescue: monotone via passive imputation
    if (k >= 2) {
      meth["R1"] <- "logreg"
    }

    if (k > 2) {
      for (j in 2:(k - 1)) {
        meth[paste0("R", j)] <-
          paste0("~ I(pmax(R", j - 1, ", R", j, "))")
      }
    }

    # age and trt not imputed
    meth[c("age", "trt")] <- ""

    # Predictor matrix
    pred <- make.predictorMatrix(dat[vars_imp])

    diag(pred) <- 0

    # trt not used as predictor (impute by arm)
    pred[, "trt"] <- 0

    # baseline always predictor
    pred[, "y0"] <- 1

    # Imputation by treatment group
    imp_list <- vector("list", 2)

    for (g in 0:1) {

      dat_g <- dat %>%
        filter(trt == g) %>%
        select(all_of(vars_imp))

      # ✅ Convert rescue variables to factors BEFORE imputation
      if (length(vars_R) > 0) {
        dat_g[vars_R] <- lapply(dat_g[vars_R],
                                function(x) factor(x, levels = c(0, 1)))
      }

      imp_list[[g + 1]] <- mice(
        dat_g,
        m = m,
        method = meth,
        predictorMatrix = pred,
        maxit = maxit,
        seed = seed + g,
        printFlag = FALSE
      )
    }

    # Recombine completed datasets
    imp_full <- lapply(1:m, function(i) {
      bind_rows(
        complete(imp_list[[1]], i),
        complete(imp_list[[2]], i)
      )
    })

    # Fit ANCOVA model in each imputed dataset
    fit_models <- lapply(imp_full, function(d) {
      d$chg <- d[[paste0("y", k)]] - d$y0
      lm(chg ~ trt + age + y0, data = d)
    })

    # Pool results
    imp_obj <- as.mira(fit_models)
    pooled <- pool(imp_obj)

    sum_pooled <- summary(
      pooled,
      conf.int = TRUE,
      conf.level = ci_level
    )

    trt_row <- sum_pooled[sum_pooled$term == "trt", ]

    list(
      coef     = trt_row$estimate,
      p        = trt_row$p.value,
      ci_lower = trt_row[[grep("%", names(trt_row))[1]]],
      ci_upper = trt_row[[grep("%", names(trt_row))[2]]]
    )
  }
}

Design <- assumptions_diabetes_rescue()
condition <- Design[1, ]

dat <- generate_diabetes_rescue(condition)

analyse_mi <- analyse_diabetes_rescue_mi()

analyse_mi(condition, dat)
