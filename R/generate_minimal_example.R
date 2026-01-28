#' Generate dataset for minimal example
#'
#' @param condition row of Design dataset
#' @param fixed_objects list of parameters that are fixed across simulations
#'
#' @returns a simulated dataset
#' @export
#'
#' @examples
#' Design <- assumptions_minimal_example() |>
#' true_summary_statistics_minimal_example()
#'
#' dat <- generate_minimal_example(Design[1,])
#' head(dat)
generate_minimal_example <- function(condition, fixed_objects=NULL){
  data.frame(
    group = c(rep(1, condition$n), rep(0, condition$n)),
    y = c(rnorm(condition$n, condition$mean1), rnorm(condition$n, condition$mean0))
  )
}

#' Calculate true summary statistics value for generate_minimal_example
#'
#' @param Design
#'
#' @returns For true_summary_statistics_minimal_example: the Design tibble with added column eff_size
#' @export
#' @describeIn generate_minimal_example calculate true effect size
true_summary_statistics_minimal_example <- function(Design){
  Design$eff_size <- Design$mean1 - Design$mean0

  Design
}

#' Create an empty assumtions tibble for generate_minimal_example
#'
#' @param print print code to generate parameter set?
#'
#' @returns For assumptions_mimimal_example: a design tibble with default values invisibly
#' @export
#' @describeIn generate_minimal_example generate default assumptions `tibble`
assumptions_minimal_example <- function(print=interactive()){
  skel <- "params_scenarios_grid(
    n = c(50, 100), # n, pre group
    mean1 = c(1,0), # mean, group 1
    mean0 = 0       # mean, group 0
  )"

  if(print){
    cat(skel)
  }

  invisible(
    skel |>
      str2expression() |>
      eval()
  )
}

