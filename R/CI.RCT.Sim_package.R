#' CI.RCT.Sim: simulate scenarios in which causal inference methods are applied in randomized controlled trials
#'
#' Package to Facilitate a Simulation Study on Causal Inference Methods in RCTs
#'
#' @keywords internal
#' @name CI.RCT.Sim
#' @title CI.RCT.Sim
"_PACKAGE"

## usethis namespace: start
#' @import SimDesign
NULL

# declaring varibles to avoid R CMD check notes.
# those are mostly column names that occur in with, within, subset functions,
# dplyr verbs and ggplot calls.
globalVariables(c(
  "trt", "R", "id", "Visit", "d", "rescue_start", "y0", "age", "exposure",
  "hba1c", "hba1c_0", "rescue", "visit", "y", "R0"
))
