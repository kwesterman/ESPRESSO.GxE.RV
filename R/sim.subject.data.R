#' 
#' @title Simulates the individual effect related to heterogeneity in baseline disease risk
#' @description The variation in baseline disease risk is assumed to be normally distributed 
#' on a logistic scale. If this parameter is set to 10, the implication is that a 'high risk' 
#' subject (someone at the upper 95 percent entile of population risk) is, all else being equal, 
#' at 10 times the offs of developing disease compared to someone else who is at 'low risk' (at 
#' the lower 5 percent centile of population risk).
#' @param num.obs number of observations to simulate.
#' @param baseline.OR baseline odds ratio for subject on 95 percent population centile versus 5 
#' percentile. This parameter reflects the heterogeneity in disease risk arising from determinantes 
#' that have not been measured or have not been included in the model.
#' @return a numerical vector.
#' @keywords internal
#' @author Gaye A.
#'
sim.subject.data <- function (num.obs=NULL, baseline.OR=NULL){
  
  numobs <- num.obs
  baseline.odds <- baseline.OR
  
  # CONVERT BASELINE ODDS RATIO FROM 5th TO 95th PERCENTILES INTO THE
  # CORRESPONDING VARIANCE FOR A NORMALLY DISTRIBUTED RANDOM EFFECT 
  baseline.variance <- (log(baseline.odds)/(2*qnorm(0.95)))^2
  
  # CREATE NORMALLY DISTRIBUTED RANDOM EFFECT VECTOR
  # WITH APPROPRIATE VARIANCE ON SCALE OF LOG-ODDS
  subject.effect <- rnorm(numobs,0,sqrt(baseline.variance))
  
  # RETURN A VECTOR 
  output <- subject.effect
}
