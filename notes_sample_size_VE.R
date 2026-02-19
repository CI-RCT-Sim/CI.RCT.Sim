

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



vaccine_scenario_set_samplesize <- function(Design){
  set_samplesize_rowwise <- function(condition){
    ns <- sample_size_formula_nauta(0.025, 1-condition$rr_ps, miniPCH::ppch(365, c(0, 14), c(0, condition$lambda_post)), 0.3, 0.9, 1)
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


Design <- vaccine_scenario(print=FALSE) |>
  vaccine_scenario_set_gamma_0() |>
  vaccine_scenario_set_true_eff() |>
  vaccine_scenario_set_samplesize()

Design
