# ------------------------------------------------------------------------------
# OBJECTIVE FUNCTION FOR MINIMIZATION OF N WITH PENALTIES FOR MISLEADING EVIDENCE (Min Mean)
# ------------------------------------------------------------------------------
# COLLAPSING BOUNDARIES (WEIBULL DISTRIBUTION)
# ------------------------------------------------------------------------------

# params: vector consisting of 
#'@param a initial value (at x = 0) for H0 boundary
#'@param asympt asymptotic value for H0 boundary
#'@param lambda scale parameter for H0 boundary
#'@param k shape parameter for H0 boundary
#'!!! params a and asympt are on the log scale and should be defined for H1
#'
#'Other arguments
#'@param sim.H0 simulated sequential data under H0
#'@param sim.H1 simulated sequential data under H1 (some ES != 0)
#'@param max.FP maximum allowed rate of false positive evidence
#'@param max.FN maximum allowed rate of false negative evidence
#'@param BF.column column number of the Bayes factors in the sim objects

#'@import dplyr

minimize.meanN.r <- function(params, sim.H0, sim.H1, max.FP = 0.05, max.FN = 0.2, BF.column = 5){
  
  a.H1 <- params[1]
  a.H0 <- -params[1]
  asympt.H1 <- params[2]
  asympt.H0 <- -params[2]
  lambda.H1 <- params[3]
  lambda.H0 <- params[3]
  k.H1 <- params[4]
  k.H0 <- params[4]
  
  require(dplyr)
  
  # Define boundaries for each N in the sim objects
  sim.H0$logb.H0 <- weibull(sim.H0$n, a = a.H0, asympt = asympt.H0, lambda = lambda.H0, k = k.H0)
  sim.H0$logb.H1 <- weibull(sim.H0$n, a = a.H1, asympt = asympt.H1, lambda = lambda.H1, k = k.H1)
  
  
  sim.H1$logb.H0 <- weibull(sim.H1$n, a = a.H0, asympt = asympt.H0, lambda = lambda.H0, k = k.H0)
  sim.H1$logb.H1 <- weibull(sim.H1$n, a = a.H1, asympt = asympt.H1, lambda = lambda.H1, k = k.H1)
  
  # Get Results under H0
  res.H0 <- sim.H0 %>%
    filter(.[[BF.column]] >= logb.H1 | .[[BF.column]] <= logb.H0) %>%
    group_by(id) %>%
    filter(row_number()==1) %>%
    ungroup()
  
  # Compute FP Evidence Rate and Mean under H0
  mean.H0 <- mean(res.H0$n)
  FP.rate <- nrow(res.H0 %>% filter(.[[BF.column]] >= logb.H1))/nrow(res.H0)
  
  # Get Results under H1
  res.H1 <- sim.H1 %>%
    filter(.[[BF.column]] >= logb.H1 | .[[BF.column]] <= logb.H0) %>%
    group_by(id) %>%
    filter(row_number()==1) %>%
    ungroup()
  
  # Compute FN Evidence Rate and Mean under H1
  mean.H1 <- mean(res.H1$n)
  FN.rate <- nrow(res.H1 %>% filter(.[[BF.column]] <= logb.H0))/nrow(res.H1)
  
  # Compute Penalty
  if (FP.rate > max.FP | FN.rate > max.FN) {
    penalty <- 10000 + (1/mean(exp(res.H0$logb.H0)) + mean(exp(res.H0$logb.H1)) +
                          1/mean(exp(res.H1$logb.H0)) + mean(exp(res.H1$logb.H1)))/2
  } else if (params[1] < params[2]){
    penalty <- 10000 + (1/mean(exp(res.H0$logb.H0)) + mean(exp(res.H0$logb.H1)) +
                          1/mean(exp(res.H1$logb.H0)) + mean(exp(res.H1$logb.H1)))/2
  } else {
    penalty <- 0
  }
  
  # Compute Result
  N <- (mean.H0 + mean.H1)/2 + penalty
  
  return(N)
  
}

