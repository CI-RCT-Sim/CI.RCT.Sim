devtools::load_all()

condition <- diabetes_scenario() |>
  diabetes_scenario_set_truevalues() |>
  _[1, ]


dat <- generate_diabetes(
  condition,
  fixed_objects=NULL
)


iteration <- 2

# investigate mi ----------------------------------------------------------

Rprof(paste0("mi_", iteration, ".Rprof"))

analyse_diabetes_mi(strategy = "hypothetical")(condition, dat, fixed_objects=NULL)

Rprof(NULL)
summaryRprof(paste0("mi_", iteration, ".Rprof"))

# investigate gcomp -------------------------------------------------------

Rprof(paste0("gcomp_", iteration, ".Rprof"))

analyse_diabetes_gcomputation()(condition, dat, fixed_objects=NULL)

Rprof(NULL)
summaryRprof(paste0("gcomp_", iteration, ".Rprof"))
