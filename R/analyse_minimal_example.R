#' Analyse dataset from minimal example scenario with linear regression
#'
#' @param ci_level the confidence level for the CIs (defaults to 0.95)
#'
#' @returns an analyse function that returns a list with the elements
#'  * `p` the p-value of the F-test
#'  * `coef` the point estimate
#'  * `ci_lower` the lower CI limit
#'  * `ci_upper` the upper CI limit
#' @export
#'
#' @importFrom stats lm confint confint.lm anova coefficients
#'
#' @examples
#' Design <- assumptions_minimal_example() |>
#'   true_summary_statistics_minimal_example()
#'
#' condition <- Design[1, ]
#'
#' dat <- generate_minimal_example(condition)
#'
#' my_analyse_lm <- analyse_minimal_example_lm(ci_level=0.9)
#' my_analyse_lm(condition, dat)
analyse_minimal_example_lm <- function(ci_level=0.95){
  function(condition, dat, fixed_objects = NULL){
    model <- lm(y~group, dat)
    CI <- confint(model, level=ci_level)

    list(
      p = anova(model)$`Pr(>F)`[1],
      coef = coefficients(model)["group"],
      ci_lower = CI["group", 1],
      ci_upper = CI["group", 2]
    )
  }
}



#' Analyse dataset from minimal example scenario with t-test
#'
#' @returns an analyse function that returns a list with the elements
#'  * `p` the p-value
#' @export
#' @importFrom stats t.test
#'
#' @examples
#' Design <- assumptions_minimal_example() |>
#'   true_summary_statistics_minimal_example()
#'
#' condition <- Design[1, ]
#'
#' dat <- generate_minimal_example(condition)
#'
#' my_analyse_t  <- analyse_minimal_example_t()
#' my_analyse_t(condition, dat)
analyse_minimal_example_t <- function(){
  function(condition, dat, fixed_objects = NULL){
    ttest <- t.test(dat$y[dat$group==1], dat$y[dat$group==0])

    list(
      p = ttest$p.value
    )
  }
}


