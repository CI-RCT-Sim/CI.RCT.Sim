## Sample size justification and validation of corresponding algorithms against gsDesign
library(tidyverse)
devtools::load_all()

#' Compute sample size requirements for VE trial using Lachin-Foulkes method
#'
#' This is useful as it provides some information on the number of expected events and total sample size required to achieve it, starting with assumptions about case incidence, dropout and follow-up
#'
#' @param VE assumed vaccine efficacy
#' @param VE0 lower bound required for demonstrating VE
#' @param r randomization ratio (1 corresponds to balanced randomization)
#' @param alpha one-sided alpha
#' @param beta 1-Power
#' @param lambdaC 1 month control arm case incidence
#' @param dropoutRate dropout rate
#' @param enrollDuration time to enroll target subject number
#' @param minfup minimum follow-up
#'
#' @return an nSurv object giving details on the trial characteristics
#' @export
#'
#' @examples
#' # Demonstrate better than 0.3 VE for a vaccine with VE of .7 with power .9
#' n_lachin_foulkes(VE = .7, VE0 = .3, lambdaC = 4/3000, dropoutRate = 0, enrollDuration = 3, minfup = 1)
n_lachin_foulkes <- function(VE=.6,VE0=.3,
                             r = 1,
                             alpha = 0.025,
                             beta = 0.1,
                             lambdaC = 0.01,
                             dropoutRate = 0.05,
                             enrollDuration = 3,
                             minfup = 3
){
  hr = 1-VE
  hr0 = 1-VE0
  gsDesign::nSurv(hr = hr, hr0 = hr0,lambdaC = lambdaC,
                  eta = dropoutRate,T = enrollDuration+minfup,minfup = minfup,alpha=alpha,beta=beta,r=r)
}

n_lachin_foulkes(VE = 0.7,
                 VE0 = .3,
                 lambdaC = 1/2000,
                 dropoutRate = 0,
                 enrollDuration = 0.1,
                 minfup = 6,
                 beta = 0.2)$n

## Example code to obtain sample size estimates
vaccine_samplesizes <- vaccine_scenario() |> vaccine_scenario_set_beta_A1_relative() |> vaccine_scenario_set_gamma_0() |> vaccine_scenario_set_samplesize()

vaccine_samplesizes |>
  rowwise() |>
  mutate(n_check = n_lachin_foulkes(VE = 1-exp(beta_A2),
                                    VE0 = .3,
                                    lambdaC = lambda_post,
                                    dropoutRate = 0,
                                    enrollDuration = 0.1,
                                    minfup = follow_up-2,
                                    beta = 0.2)$n,
         n_sim = n_trt+n_ctrl,
         diff = round((n_sim-n_check)/n_check*100,2)) -> out

## To check vaccine efficacies use:
out |> mutate(across(starts_with('beta'),~ exp(.)),
              `1/monthly_incidence` = 1/(1-exp(-lambda_post*(30/7)))) |> select(`1/monthly_incidence`,beta_A1,beta_A2,n_sim,n_check,diff)

