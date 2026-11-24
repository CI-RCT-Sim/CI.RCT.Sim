#' Analyse data set with G-computation
#'
#' @return A function that, when called with `condition` and `dat`, returns a list with:
#' * `coef`      estimated difference in mean change in HbA1c between treatment groups
#' * `se`        standard error for coef
#' * `p`         p-value for coef
#' * `ci_lower`  lower bound of 95% confidence interval for coef
#' * `ci_upper`  upper bound of 95% confidence interval for coef
#'
#' @export
#'
#' @importFrom dplyr mutate arrange group_by ungroup lead lag filter select
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr `%>%`
#'
#' @details
#' This function implements G-computation to estimate the
#' average treatment effect (ATE) on the **change in HbA1c** under a hypothetical scenario
#' where rescue medication is not available, regardless of treatment discontinuation.
#'
#' Two scenarios are considered depending on `setup`:
#' - `setup = 1`: the effect of rescue goes on top of the effect of the assigned treatment
#' - `setup = 0`: once rescue starts, treatment is deterministically recoded to control
#'
#' The estimand is:
#' - **Endpoint**: Change in HbA1c at final visit k
#' - **Treatments**: Experimental vs. control
#' - **Summary measure**: Difference in means
#' - **Population**: People with type 2 diabetes
#' - **Intercurrent events**: Use of rescue medication and treatment discontinuation
#' - **Intercurrent event strategy**:
#'   - *Hypothetical strategy*: Assume rescue medication was not available (i.e., rescue = 0 always)
#'   - *Treatment policy strategy*: Ignore treatment discontinuation (i.e., continue treatment in counterfactual)
#'
# Data are reshaped to long format. Baseline HbA1c, change from baseline,
# lagged HbA1c, and the final outcome y_k are constructed, and the final visit
# row is removed so that visit k - 1 contains the full covariate history needed
# to model the end-of-follow-up outcome.
#
# The fitted models are:
#   hba1c ~ trt + R + lag1_hba1c + age
#   y_k   ~ trt + R + lag1_hba1c + age + hba1c_0
#
# Counterfactual trajectories are then simulated separately for each treatment arm,
# starting from the observed baseline covariate distribution and imposing R = 0
# at all post-baseline visits. HbA1c is simulated forward under the intervention,
# and the final change in HbA1c is simulated from the outcome model.
# The treatment effect is the difference in mean simulated final change between arms.
#
# A separate rescue model is not used in the simulation step because the intervention
# is deterministic: rescue is fixed to 0 at all times. Thus, for this specific
# no-rescue estimand, the post-intervention distribution can be simulated directly
# from the HbA1c and final outcome models evaluated at R = 0. A rescue model would
# be needed for a more general g-formula implementation of the natural course or
# alternative rescue strategies.
#
# Uncertainty is quantified by bootstrapping (200 samples, percentile method) the entire procedure:
# resample subjects, rebuild the long data, refit the models, resimulate the
# no-rescue counterfactual outcomes, and recompute the treatment contrast.
#'
#' The null hypothesis of no treatment effect is rejected if the 95% CI does not include 0.
#'
#' @examples
#' \donttest{
#' Design <- diabetes_scenario()[1, ] |>
#'   diabetes_scenario_set_truevalues()
#'
#' dat <- generate_diabetes(Design)
#'
#' analyse_diabetes_gcomputation()(Design, dat)
#' }
analyse_diabetes_gcomputation <- function() {
  function(condition, dat, fixed_objects = NULL,B=200) {
    data<-dat
    k<-condition$k
    estim_fun<-function(data) {
      mean_age<-mean(data$age)
      mean_y<-mean(data$y0)

      group<-c(0,1)
      yfinal<-rep(NA,2)
      for(u in 1:2) {
        yhat<-rep(NA,k)

        D<-data[data$trt==group[u],]

        X<-matrix(c(1,mean_y,0,mean_age),nrow=1)
        i<-1
        frml<-formula(paste("y",i,"~ y",i-1,"+R",i-1,"+age",sep=""))
        mod<-lm(frml,data=D)
        beta<-coef(mod)
        beta[is.na(beta)]<-0 #because if there is no rescue, the coefficient is NA
        yhat[1]<-X%*%beta

        for(i in 2:k) {
          X<-matrix(c(1,yhat[i-1],0,mean_y,mean_age),nrow=1)
          frml<-formula(paste("y",i,"~ y",i-1,"+R",i-1,"+y0+age",sep=""))
          mod<-lm(frml,data=D)
          beta<-coef(mod)
          beta[is.na(beta)]<-0 #because if there is no rescue, the coefficient is NA
          yhat[i]<-X%*%beta
        }
        yfinal[u]<-yhat[k]
      }
      est<-yfinal[2]-yfinal[1]
      est
    }

    est<-estim_fun(data)

    #Bootsratp
    #B<-100
    ids<-1:dim(data)[1]
    est_b<-rep(NA,B)
    for(j in 1:B) {
      ind<-sample(ids,replace=TRUE)
      data_b<-data[ind,]
      est_b[j]<-estim_fun(data_b)
    }
    SE<-sd(est_b,na.rm=TRUE)
    Q<-qt(0.975,df=B-1)
    KI<-est+c(-1,1)*Q*SE
    KI
    z<-est/SE
    p<-pnorm(z)
    p
    list(
      coef = est,
      se = SE,
      p = p,
      ci_lower = KI[1],
      ci_upper = KI[2]
    )
  }
}

