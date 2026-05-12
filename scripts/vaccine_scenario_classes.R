library(SimNPH)
library(tidyverse)

vaccine_scenario_defaults()

scenario_A1 <- vaccine_scenario_base_defaults()
scenario_A1$lambda_post <- vaccine_scenario_defaults()$lambda_post
scenario_A1$effect_before_d2 <- TRUE
scenario_A1$beta_A2 <- vaccine_scenario_defaults()$beta_A2

vaccine_scenario_A1 <- vaccine_scenario(scenario_defaults = \() scenario_A1)

scenario_A2 <- vaccine_scenario_base_defaults()
scenario_A2$effect_before_d2 <- TRUE
scenario_A2$beta_A2 <- vaccine_scenario_defaults()$beta_A2
scenario_A2$gamma_A <- vaccine_scenario_defaults()$gamma_A

vaccine_scenario_A2 <- vaccine_scenario(scenario_defaults = \() scenario_A2)  # maybe filter gamma_A != 0

scenario_B1 <- vaccine_scenario_base_defaults()
scenario_B1$effect_before_d2 <- TRUE
scenario_B1$p_V <- vaccine_scenario_defaults()$p_V[3]
scenario_B1$beta_A2 <- vaccine_scenario_defaults()$beta_A2
scenario_B1$gamma_A <- vaccine_scenario_defaults()$gamma_A
scenario_B1$gamma_V <- vaccine_scenario_defaults()$gamma_V
scenario_B1$beta_V <- vaccine_scenario_defaults()$beta_V

vaccine_scenario_B1 <- vaccine_scenario(scenario_defaults = \() scenario_B1)
## maybe B2 with p_V = 0.1

scenario_C1 <- vaccine_scenario_base_defaults()
scenario_C1$effect_before_d2 <- TRUE
scenario_C1$p_W <- vaccine_scenario_defaults()$p_W[3]
scenario_C1$beta_A2 <- vaccine_scenario_defaults()$beta_A2
scenario_C1$gamma_A <- vaccine_scenario_defaults()$gamma_A
scenario_C1$gamma_W <- vaccine_scenario_defaults()$gamma_W
scenario_C1$gamma_AW <- vaccine_scenario_defaults()$gamma_AW
scenario_C1$beta_W <- vaccine_scenario_defaults()$beta_W
scenario_C1$beta_AW <- vaccine_scenario_defaults()$beta_AW

vaccine_scenario_C1 <- vaccine_scenario(scenario_defaults = \() scenario_C1)
## mabe C2 with p_V = 0.1 or a bit more variation in the combination of gamma_AW and beta_AW

scenario_D1 <- vaccine_scenario_base_defaults()
scenario_D1$effect_before_d2 <- TRUE
scenario_D1$p_W <- vaccine_scenario_defaults()$p_W[3]
scenario_D1$beta_A2 <- vaccine_scenario_defaults()$beta_A2
scenario_D1$gamma_A <- vaccine_scenario_defaults()$gamma_A
scenario_D1$gamma_W <- vaccine_scenario_defaults()$gamma_W
scenario_D1$gamma_AW <- vaccine_scenario_defaults()$gamma_AW
scenario_D1$beta_W <- vaccine_scenario_defaults()$beta_W
scenario_D1$beta_AW <- vaccine_scenario_defaults()$beta_AW
scenario_D1$p_V <- vaccine_scenario_defaults()$p_V[3]
scenario_D1$gamma_V <- vaccine_scenario_defaults()$gamma_V
scenario_D1$beta_V <- vaccine_scenario_defaults()$beta_V

vaccine_scenario_D1 <- vaccine_scenario(scenario_defaults = \() scenario_D1)
## likely D1 with more combinations of factors could be interesting


## overall compliance to complicated scenarios
scenario_extra$overall_compliance <- vaccine_scenario_defaults()$overall_compliance
## effect before D2 false to complicated scenarios
scenario_extra$effect_before_d2 <- vaccine_scenario_defaults()$effect_before_d2

scenario_extra <- vaccine_scenario_tweak(
  overall_compliance = c(0.95,0.5),
  effect_before_d2 = c(T,F),
  gamma_A = vaccine_scenario_defaults()$gamma_A,
  gamma_W = c(-0.8,0.8),
  gamma_V = c(0.5,-0.5),
  gamma_AW = c(-0.3,0.3),
  p_V = 0.3,
  p_W = 0.3,
  beta_V = vaccine_scenario_defaults()$beta_V[1]*c(1,-1),
  beta_W = vaccine_scenario_defaults()$beta_W[1]*c(1,-1),
  beta_AW = vaccine_scenario_defaults()$beta_AW[1]
)
scenario_extra

# Compliance:             C[+,-] + better compliance, - worse compliance
# Prognosis:              P[+,-] + less infection risk, - more infection risk
# Treatment modification: T[+,-] + better efficacy, - worse efficacy
names_extra <- c("CA-.CV+.CW-.CAW-.PV-.PW-.TA+.TAW+",
                 "CA-.CV+.CW-.CAW-.PV-.PW-.TA+.TAW+.OC50",
                 "CA-.CV+.CW+.CAW-.PV-.PW-.TA+.TAW+",
                 "CA-.CV-.CW-.CAW-.PV-.PW-.TA+.TAW+",
                 "CA0.CV+.CW-.CAW-.PV-.PW-.TA+.TAW+",
                 "CA-.CV+.CW-.CAW+.PV-.PW-.TA+.TAW+",
                 "CA-.CV+.CW-.CAW-.PV+.PW-.TA+.TAW+",
                 "CA-.CV+.CW-.CAW-.PV-.PW+.TA+.TAW+",
                 "CA-.CV+.CW-.CAW-.PV-.PW-.TA+.TAW+.D1F",
                 )
