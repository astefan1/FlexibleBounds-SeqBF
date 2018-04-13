# ------------------------------------------------------------------------------
# OBJECTIVE FUNCTION FOR MINIMIZATION OF N WITH PENALTIES FOR MISLEADING EVIDENCE
# ------------------------------------------------------------------------------
# LINEAR NON-SYMMETRIC BOUNDARIES, CRITERION: MEANS
# ------------------------------------------------------------------------------

#'@param bounds vector: c(lower boundary, upper boundary)
#'@param sim.H0 simulated sequential data under H0
#'@param sim.H1 simulated sequential data under H1 (some ES != 0)
#'@param max.FP maximum allowed rate of false positive evidence
#'@param max.FN maximum allowed rate of false negative evidence
#'@param BF.column column number of the Bayes factors in the sim objects

#'@import dplyr

minimize.meanN.linear <- function(bounds, sim.H0, sim.H1, max.FP = 0.05, max.FN = 0.2, BF.column = 5){
  
  require(dplyr)
  
  # logarithmize boundaries
  log.a <- log(bounds[2]) #upper
  log.b <- log(bounds[1]) #lower
  
  # Get Results under H0
  res.H0 <- sim.H0 %>%
    filter(.[[BF.column]] >= log.a | .[[BF.column]] <= log.b) %>%
    group_by(id) %>%
    filter(row_number()==1) %>%
    ungroup()
  
  # Compute FP Evidence Rate and Mean under H0
  mean.H0 <- mean(res.H0$n)
  FP.rate <- nrow(res.H0 %>% filter(.[[BF.column]] >= log.a))/nrow(res.H0)
  
  # Get Results under H1
  res.H1 <- sim.H1 %>%
    filter(.[[BF.column]] >= log.a | .[[BF.column]] <= log.b) %>%
    group_by(id) %>%
    filter(row_number()==1) %>%
    ungroup()
  
  # Compute FN Evidence Rate and Mean under H1
  mean.H1 <- mean(res.H1$n)
  FN.rate <- nrow(res.H1 %>% filter(.[[BF.column]] <= log.b))/nrow(res.H1)
  
  # Compute Penalty
  if (FP.rate > max.FP | FN.rate > max.FN) {
    penalty <- 10000 + 1/bounds[1] + bounds[2]
  } else {
    penalty <- 0
  }
  
  # Compute Result
  N <- (mean.H0 + mean.H1)/2 + penalty
  
  return(N)
  
}