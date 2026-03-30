#' Generate Data for Vaccine Scenario
#' @param print print code to generate parameter set?
#'
#' @return For `vaccine_scenario`: a design tibble with default values invisibly
#'
#' @details vaccine_scenario generates a default design `data.frame`
#'   for use with `generate_vaccine`. If print is `TRUE` code to produce
#'   the template is also printed for copying, pasting and editing by the user.
#'   (This is the default when run in an interactive session.)
#' `vaccine_scenario_grid()` constructs a simulation design by evaluating the
#' named parameter list returned by `scenario_defaults` via
#' `params_scenarios_grid()`. It then derives the parameter `beta_A1` from
#' `beta_A2` and `effect_before_d2`.
#'
#' If `print = TRUE`, the function also prints code that reproduces the same
#' parameter grid specification. This is intended to make it easy to copy,
#' modify, and reuse scenario settings.
#'
#' The default settings generate a factorial grid across the vector-valued
#' elements returned by `vaccine_scenario_defaults()`. A custom defaults function
#' can be supplied to generate an alternative grid.
#'
#' @export
#' @describeIn vaccine_scenario Generate default assumptions `data.frame`
#'
#' @examples
#' Design <- vaccine_scenario() |>
#'   vaccine_scenario_set_beta_A1_relative() |>
#'   vaccine_scenario_set_gamma_0() |>
#'   vaccine_scenario_set_true_eff() |>
#'   vaccine_scenario_set_samplesize()
#' Design
#'
#' generate_vaccine(Design[1,])
vaccine_scenario <- function(print = interactive(),scenario_defaults = vaccine_scenario_defaults) {

  defaults <- scenario_defaults()

  if (print) {
    code <- paste0(
      "params_scenarios_grid(\n",
      paste(
        sprintf("  %s = %s", names(defaults),
                vapply(defaults, function(x) paste(deparse(x), collapse = " "), character(1))
        ),
        collapse = ",\n"
      ),
      "\n)"
    )
    cat(code)
  }

  out <- do.call(params_scenarios_grid, defaults)

  invisible(out)
}

#' @param condition Row of parameters dataset
#' @param fixed_objects other objects passed to simulation runs
#'
#' @describeIn vaccine_scenario Simulate a dataset from the vaccine scenario
#' @importFrom stats binomial runif rbinom
#' @importFrom miniPCH rpch
#'
#' @returns For `generate_vaccine`: a simulated dataset
#' @export
generate_vaccine <- function(condition, fixed_objects = list(include_unobserved=FALSE)){
  N <- condition$n_ctrl + condition$n_trt
  V <- rbinom(N, 1, condition$p_V)
  W <- rbinom(N, 1, condition$p_W)
  A <- c(rep(1, condition$n_trt), rep(0, condition$n_ctrl))

  # regression model for compliance propensity
  p <- function(A_){
    binomial()$linkinv(
      condition$gamma_0  +
        condition$gamma_W * W +
        condition$gamma_V * V +
        condition$gamma_A * A_ +
        condition$gamma_AW * A_ * W
    )
  }

  # latent variable
  U <- runif(N)
  # counterfactual compliances
  C_1 <- as.numeric(U < p(1))
  C_0 <- as.numeric(U < p(0))
  # observed compliance
  C <- as.numeric(U < p(A))

  # beta_AW^(1) == beta_AW^(2) == beta_AW according to table 4.1. in the protocol
  theta_1 <- condition$beta_A1 + condition$beta_AW * W
  theta_2 <- condition$beta_A2 + condition$beta_AW * W
  theta_early <- theta_1
  theta_late  <- theta_1 + C * (theta_2 - theta_1)

  t_ <- c(0, condition$dose_interval)
  lambda_0 <- diag(c(0, condition$lambda_post))

  # under potential treatment allocation to control
  # time to event hazard
  lambda_a0 <- cbind(
    exp(condition$beta_V * V + condition$beta_W * W),
    exp(condition$beta_V * V + condition$beta_W * W)
  ) %*% lambda_0
  # time to event outcome
  T_a0 <- sapply(1:N, \(i){
    rpch(1, t=t_, lambda = lambda_a0[i,])
  })
  # binary outcome
  Y_a0 <- T_a0 < condition$follow_up
  T_a0 <- pmin(T_a0, condition$follow_up)

  # under potential treatment allocation to active
  # time to event hazard
  lambda_a1 <- cbind(
    exp(condition$beta_V * V + condition$beta_W * W + 1 * theta_early),
    exp(condition$beta_V * V + condition$beta_W * W + 1 * theta_late)
  ) %*% lambda_0
  # time to event outcome
  T_a1 <- sapply(1:N, \(i){
    rpch(1, t=t_, lambda = lambda_a1[i,])
  })
  # binary outcome
  Y_a1 <- T_a1 < condition$follow_up
  T_a1 <- pmin(T_a1, condition$follow_up)


  # time to event outcome
  T_ <- ifelse(A==0, T_a0, T_a1)
  # binary outcome
  Y <- T_ < condition$follow_up
  # time to event
  T_ <- pmin(T_, condition$follow_up)

  if(fixed_objects$include_unobserved){
    data.frame(
      V=V,
      W=W,
      trt=A,
      C=C,
      t=T_,
      evt=Y,
      p=p(A),
      C_1=C_1,
      C_0=C_0,
      evt_1=Y_a1,
      evt_0=Y_a0,
      t_1=T_a1,
      t_0=T_a0
    )
  } else {
    data.frame(
      V=V,
      W=W,
      trt=A,
      C=C,
      t=T_,
      evt=Y
    )
  }
}



#' @describeIn vaccine_scenario Default parameter values for the full vaccine simulation scenario grid
#'
#' @return For `vaccine_scenario_defaults:` A named list of default parameter
#' values for use with `vaccine_scenario_grid()`.
#'
#' @details
#' This function returns the default parameter settings used to generate the
#' standard simulation design for `generate_vaccine`. In contrast to
#' `vaccine_scenario_base_defaults()`, several entries are vectors and therefore
#' define a multi-row parameter grid when passed to `params_scenarios_grid()`.
#'
#' The returned list is primarily intended as an internal or user-modifiable
#' source of defaults that can be supplied via the `scenario_defaults` argument
#' of `vaccine_scenario_grid()`.
#'
#' @examples
#' defaults <- vaccine_scenario_defaults()
#' str(defaults)
#' defaults$beta_A2
#'
#' @export
vaccine_scenario_defaults <- function() {
  list(
    p_V = c(0, 0.1, 0.3), # probability for binary covariate prognostic for ICE and infection risk
    p_W = c(0, 0.1, 0.3), # probability for binary covariate prognostic for ICE and infection risk and modifying treatment effect
    lambda_post = -log(1 - (1 / c(1000, 500, 2000))) / (30 / 7), # weekly baseline infection hazard after 14 days
    overall_compliance = c(0.95, 0.5), # used to calibrate gamma0
    gamma_W = c(-0.8, 0), # compliance modifier predictive covariate
    gamma_V = c(0.5, 0), # compliance modifier prognostic covariate
    gamma_A = c(-0.357, 0), # compliance modifier treatment group
    gamma_AW = c(-0.3, 0), # compliance modifier treatment interaction predictive covariate
    beta_V = c(log(1.5), 0), # prognostic effect prognostic covariate
    beta_W = c(log(1.2), 0), # prognostic effect predictive covariate
    effect_before_d2 = c(TRUE, FALSE), # indicator whether one dose has an effect
    beta_A2 = log(1 - c(0.7, 0, 0.3, 0.5, 0.9)), # vaccine efficacy with full vaccination
    beta_AW = c(log(0.8), 0), # effect modifier predictive covariate
    follow_up = 26, # weeks
    dose_interval = 2 # weeks
  )
}

#' @describeIn vaccine_scenario Default parameter values for a single baseline vaccine simulation scenario
#'
#' @return For `vaccine_scenario_base_defaults:`
#' A named list of baseline parameter values for use with
#' `vaccine_scenario_tweak()`. Each element is a scalar and together they define
#' a single reference scenario.
#'
#' @details
#' `vaccine_scenario_base_defaults` returns a one-value-per-parameter baseline configuration for
#' `generate_vaccine`. It is intended for targeted testing and debugging, where a
#' user starts from a single reference scenario and modifies only selected
#' parameters.
#'
#' In contrast to vaccine_scenario_defaults(), the returned list does not
#' define a broad factorial grid unless the user explicitly replaces one or more
#' scalar values by vectors in vaccine_scenario_tweak().
#'
#' @examples
#' base_defaults <- vaccine_scenario_base_defaults()
#' str(base_defaults)
#' base_defaults$overall_compliance
#'
#' @export
vaccine_scenario_base_defaults <- function() {
  list(
    p_V = 0,
    p_W = 0,
    lambda_post = -log(1 - (1 / 1000)) / (30 / 7),
    overall_compliance = 0.95,
    gamma_W = 0,
    gamma_V = 0,
    gamma_A = 0,
    gamma_AW = 0,
    beta_V = 0,
    beta_W = 0,
    effect_before_d2 = FALSE,
    beta_A2 = log(1 - 0.7),
    beta_AW = 0,
    follow_up = 26,
    dose_interval = 2
  )
}


#' @describeIn vaccine_scenario Modify selected vaccine scenario parameters starting from a base scenario
#'
#' @param ... Named parameter overrides. Each name must correspond to an element
#'   returned by `scenario_defaults`. Supplied values replace the defaults. These
#'   may be scalars or vectors; vector inputs generate a small parameter grid.
#' @param scenario_defaults A function returning a named `list` of parameter
#'   defaults. By default `vaccine_scenario_base_defaults()` is used.
#'
#' @return For `vaccine_scenario_tweak`:
#' A `data.frame` containing the resulting scenario design for use with
#' `generate_vaccine`.
#'
#' @details
#' `vaccine_scenario_tweak()` is intended for targeted testing. It starts from a
#' baseline scenario returned by `scenario_defaults`, replaces only the
#' parameters supplied in `...`, and then evaluates the resulting settings with
#' `params_scenarios_grid()`.
#'
#' By default the function starts from the scalar baseline returned by
#' `vaccine_scenario_base_defaults()`, so changing only a few arguments typically
#' produces a single scenario row. If one or more supplied overrides are vectors,
#' a corresponding parameter grid is generated.
#'
#' The function validates that all supplied argument names are known parameters.
#' An error is thrown if unknown names are passed.
#'
#' As in `vaccine_scenario_grid()`, `beta_A1` is derived from `beta_A2` and
#' `effect_before_d2`.
#'
#' @examples
#' # single modified test scenario
#' Design <- vaccine_scenario_tweak(
#'   beta_A2 = log(1 - 0.5),
#'   overall_compliance = 0.8
#' )
#' Design
#'
#' # small grid around the baseline scenario
#' Design_grid <- vaccine_scenario_tweak(
#'   beta_A2 = log(1 - c(0.3, 0.5, 0.7)),
#'   effect_before_d2 = c(TRUE, FALSE)
#' )
#' Design_grid
#'
#' @export
vaccine_scenario_tweak <- function(...,scenario_defaults = vaccine_scenario_base_defaults) {
  defaults <- scenario_defaults()
  updates <- list(...)

  unknown <- setdiff(names(updates), names(defaults))
  if (length(unknown) > 0) {
    stop(
      "Unknown parameter(s): ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }

  defaults[names(updates)] <- updates

  out <- do.call(params_scenarios_grid, defaults)

  out
}

#' @param diff Difference in vaccine efficacy between dose 1 and dose 2
#'
#' @describeIn vaccine_scenario Calculate VE after dose 1 as fraction of VE after dose 2
#' @returns For `vaccine_scenario_set_beta_A1_relative`: The Design dataset with added column `beta_A1`
#' @export
vaccine_scenario_set_beta_A1_relative <- function(Design, diff=-0.2){
  Design |>
    transform(
      beta_A1 = (ifelse(
        beta_A2 == 0,
        0,
        log(1-abs(1-(exp(beta_A2) - diff))))
        * effect_before_d2)
    )
}



#' @param Design Design dataset as returned by `vaccine_scenario`
#' @param r for setting gamm_0 and sample size calculation: allocation ratio
#'
#' @describeIn vaccine_scenario Calculate gamma_0 from other parameters
#'
#' @returns For `vaccine_scenario_set_gamma_0`: The Design dataset with the added column gamma_0
#' @export
vaccine_scenario_set_gamma_0 <- function(Design, r=1){
  set_gamma_0_rowwise <- function(condition){
    # plug-in mean for W, V
    A0 <- condition$gamma_W * condition$p_W +
      condition$gamma_V * condition$p_V

    A1 <- condition$gamma_W * condition$p_W +
      condition$gamma_V * condition$p_V +
      condition$gamma_A +
      condition$gamma_AW * condition$p_W

    w1 <- 1/(1+r)
    w0 <- 1-1/(1+r)
    p <- condition$overall_compliance

    # derived with sympy
    log_solution_1 <-
      (-p*exp(A0) - p*exp(A1) + w0*exp(A0) + w1*exp(A1) - sqrt(p^2*exp(2*A0) + p^2*exp(2*A1) - 2*p^2*exp(A0 + A1) - 2*p*w0*exp(2*A0) + 2*p*w0*exp(A0 + A1) - 2*p*w1*exp(2*A1) + 2*p*w1*exp(A0 + A1) + w0^2*exp(2*A0) + 2*w0*w1*exp(A0 + A1) + w1^2*exp(2*A1)))*exp(-A0 - A1)/(p - w0 - w1)
    log_solution_2 <-
      (-p*exp(A0) - p*exp(A1) + w0*exp(A0) + w1*exp(A1) + sqrt(p^2*exp(2*A0) + p^2*exp(2*A1) - 2*p^2*exp(A0 + A1) - 2*p*w0*exp(2*A0) + 2*p*w0*exp(A0 + A1) - 2*p*w1*exp(2*A1) + 2*p*w1*exp(A0 + A1) + w0^2*exp(2*A0) + 2*w0*w1*exp(A0 + A1) + w1^2*exp(2*A1)))*exp(-A0 - A1)/(p - w0 - w1)

    # check which solution works (log of something > 0)
    if((log_solution_1 > 0) && (log_solution_2 <= 0)){
      gamma_0 <- log(log_solution_1) - log(2)
    } else if((log_solution_2 > 0) && (log_solution_1 <= 0)){
      gamma_0 <- log(log_solution_2) - log(2)
    } else {
      stop("No valid solution")
    }

    gamma_0
  }

  Design$gamma_0 <- Design |>
    split(1:nrow(Design)) |>
    sapply(set_gamma_0_rowwise)

  Design
}



#' @describeIn vaccine_scenario Calculate relative risk in the principal stratum of always compliers
#'
#' @returns `vaccine_scenario_set_true_eff`: the Design dataset with the added column rr_ps
#' @export
vaccine_scenario_set_true_eff <- function(Design){
  set_true_eff_rowwise <- function(condition){
    p <- \(a,v,w){
      binomial()$linkinv(
        condition$gamma_0 + condition$gamma_W*w + condition$gamma_V*v + condition$gamma_A*a + condition$gamma_AW*a*w
      )
    }

    m <- \(v,w){
      pmin(p(0,v,w), p(1,v,w))
    }

    Pr_vw <- \(v,w){
      pr_V <- \(v){
        ((v*condition$p_V)+(1-v)*(1-condition$p_V))
      }

      pr_W <- \(w){
        (w*condition$p_W)+(1-w)*(1-condition$p_W)
      }

      (pr_V(v) * pr_W(w) * m(v,w)) /
        sum(outer(0:1, 0:1, \(v_, w_){pr_V(v_)*pr_W(w_)*m(v_,w_)}))
    }

    pi <- Vectorize(\(y,v,w){
      t_ <- c(0, condition$dose_interval)
      lambda_0 <- diag(c(0, condition$lambda_post))
      theta_1 <- condition$beta_A1 + condition$beta_AW * w
      theta_2 <- condition$beta_A2 + condition$beta_AW * w
      theta_early <- theta_1
      theta_late  <- theta_1 + (theta_2 - theta_1)

      lambda_vw <- cbind(
        exp(condition$beta_V * v + condition$beta_W * w + y * theta_early),
        exp(condition$beta_V * v + condition$beta_W * w + y * theta_late)
      ) %*% lambda_0

      miniPCH::ppch(condition$follow_up, t_, lambda_vw)
    })

    Pr <- \(a){
      sum(outer(0:1, 0:1, \(v,w){
        pi(a,v,w) * Pr_vw(v,w)
      }))
    }

    rr_ps <- Pr(1)/Pr(0)
    rr_ps
  }

  Design$rr_ps <- Design |>
    split(1:nrow(Design)) |>
    sapply(set_true_eff_rowwise)

  Design
}





#' @describeIn vaccine_scenario Calculate sample size as if the true effect was correctly guessed
#'
#' @param alpha for sample size calculation: one-sided alpha level
#' @param CSE for sample size calculation: super-superiority margin on the VE scale
#' @param power for sample size calculation: target power
#' @param VE_H1 vaccine efficacy assumed for sample size planning
#'
#' @returns For `vaccine_scenario_set_samplesize`: the Design dataset with the added columns n_trt and n_ctrl
#' @export
vaccine_scenario_set_samplesize <- function(Design, alpha=0.025, CSE=0.3, VE_H1=0.7, power=0.8, r=1){
  sample_size_formula_nauta <- function(alpha, VE, AR0, CSE, power, r){
    theta0 <- 1-CSE
    p0 <- AR0
    p1 <- p0*(1-VE)
    Zalpha <- qnorm(1-alpha)
    Zbeta <- qnorm(power)
    A <- 1+1/r
    B <- -(theta0*(1+p0/r)+1/r+p1)
    C <- theta0*(p1+p0/r)
    R1 <- (-B-sqrt(B*B-4*A*C))/(2*A)
    R0 <- R1/theta0
    n1 <- ((
      Zalpha*sqrt(R1*(1-R1) + ((theta0^2)*r)*R0*(1-R0)) +
        Zbeta*sqrt(p1*(1-p1) + ((theta0^2)*r)*p0*(1-p0))
    ) / (p1-theta0*p0))^2
    n0 <- floor(n1/r)+1
    n1 <- floor(n1)+1

    c(n0, n1)
  }


  set_samplesize_rowwise <- function(condition){
    ns <- sample_size_formula_nauta(alpha, VE_H1, miniPCH::ppch(condition$follow_up, c(0, condition$dose_interval), c(0, condition$lambda_post)), CSE, power, r)
    condition$n_trt <- ns[2]
    condition$n_ctrl <- ns[1]
    condition
  }

  Design <- Design |>
    split(1:nrow(Design)) |>
    lapply(set_samplesize_rowwise) |>
    do.call(rbind, args=_)

  Design
}
