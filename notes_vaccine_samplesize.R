# renv::deactivate()

library(parallel)
devtools::load_all()
library(gsDesign)

cl <- makeCluster(parallel::detectCores()-1)

clusterEvalQ(cl, {
  devtools::load_all()
  library(gsDesign)
})

Design <- vaccine_scenario() |>
  vaccine_scenario_set_gamma_0() |>
  vaccine_scenario_set_true_eff() |>
  vaccine_scenario_set_samplesize()

Design |> subset(select=c(lambda_post, beta_A2, beta_A1, rr_ps, n_trt, n_ctrl))

Analyse <- function(condition, dat, fixed_objects) {
  counts <- dat |>
    subset(C==1) |>
    by(~trt, \(rows){
      data.frame(
        trt = rows$trt[1],
        evt = sum(rows$evt),
        n   = nrow(rows)
      )
    })


  ret <- data.frame(
    p = pnorm(testBinomial(counts[["1"]]$evt, counts[["0"]]$evt, counts[["1"]]$n, counts[["0"]]$n, delta=-0.3, scale="rr"))
  )
  ret
}

Summarise <- function(condition, results, fixed_objects) {
  ret <- data.frame(
    rejection = mean(results$p < 0.025)
  )
  ret
}

#-------------------------------------------------------------------



res <- runSimulation(
  design=Design,
  replications=1000,
  generate=generate_vaccine,
  analyse=Analyse,
  summarise=Summarise,
  fixed_objects = list(include_unobserved=FALSE),
  cl=cl
)

res

res$rejection

stopCluster(cl)
# renv::activate()

