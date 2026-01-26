
#' Generate Grid of Parameter Values
#'
#' @param ... vectors of parameter values
#'
#' @returns
#' Returns a data.frame of parameter values for the simulations.
#'
#' @details
#' The first value of each argument is taken to be the reference value. In
#' addition to the reference scenario the output contains one line per
#' additional value in each of the input arugments. In each of those lines only
#' one parameter is varied. Each line differs from the first line in exactly one
#' column.
#'
#' @export
#'
#' @examples
#' params_scenarios_grid(x=1:3, y=11:12, z=21)
params_scenarios_grid <- function(...){
  params_args <- list(...)
  params_ref <- purrr::map(params_args, \(x){
    x[1]
  }) |>
    as.data.frame()

  params_other <- purrr::imap(params_args, \(x, i){
    purrr::map(x[-1], \(y){
      tmp <- params_ref
      tmp[,i] <- y
      tmp
    }) |>
      purrr::list_rbind()
  }) |>
    purrr::list_rbind()

  rbind(params_ref, params_other)
}
