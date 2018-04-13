# ------------------------------------------------------------------------------
# ---------------------------- OF EVALUATION -----------------------------------
# ------------------------------------------------------------------------------

#' @import dplyr

OF.eval <- function(opt.result, sim.H0, sim.H1, method, BF.column = 5, ...){
  
  require(dplyr)
  
  match.arg(method, c("symmetric", "linear", "collapsing", "collapsing.sym"))
  
  if(method ==  "symmetric" | method == "linear"){
    
    if(method == "symmetric"){
      # logarithmize boundaries
      log.a <- log(opt.result) #upper
      log.b <- -log.a #lower 
    }
    
    if(method == "linear"){
      # logarithmize boundaries
      log.a <- log(opt.result[2]) #upper
      log.b <- log(opt.result[1]) #lower
    }
    
    # Get Results under H0 and H1
    res.H0 <- sim.H0 %>%
      filter(.[[BF.column]] >= log.a | .[[BF.column]] <= log.b) %>%
      group_by(id) %>%
      filter(row_number()==1) %>%
      ungroup()
    
    res.H1 <- sim.H1 %>%
      filter(.[[BF.column]] >= log.a | .[[BF.column]] <= log.b) %>%
      group_by(id) %>%
      filter(row_number()==1) %>%
      ungroup()
    
  } else {
    
    weibull <- function(t, a, asympt, lambda, k){
      y <- a - (1 - exp(1)**(-(t/lambda)**k))*(a-asympt)
      return(y)
    }
    
    if(method == "collapsing"){
      
      # Parameters of the weibull function
      a.H0 <- opt.result[1]
      a.H1 <- opt.result[2]
      asympt.H0 <- opt.result[3]
      asympt.H1 <- opt.result[4]
      lambda.H0 <- opt.result[5]
      lambda.H1 <- opt.result[6]
      k.H0 <- opt.result[7]
      k.H1 <- opt.result[8]
      
      # Define boundaries for each N in the sim objects
      sim.H0$logb.H0 <- log(weibull(sim.H0$n, a = a.H0, asympt = asympt.H0, lambda = lambda.H0, k = k.H0))
      sim.H0$logb.H1 <- log(weibull(sim.H0$n, a = a.H1, asympt = asympt.H1, lambda = lambda.H1, k = k.H1))
      
      
      sim.H1$logb.H0 <- log(weibull(sim.H1$n, a = a.H0, asympt = asympt.H0, lambda = lambda.H0, k = k.H0))
      sim.H1$logb.H1 <- log(weibull(sim.H1$n, a = a.H1, asympt = asympt.H1, lambda = lambda.H1, k = k.H1))
      
    }
    
    if(method == "collapsing.sym"){
      
      # Parameters of the weibull function
      a.H0 <- -opt.result[1]
      a.H1 <- opt.result[1]
      asympt.H0 <- -opt.result[2]
      asympt.H1 <- opt.result[2]
      lambda <- opt.result[3]
      k <- opt.result[4]
      
      # Define boundaries for each N in the sim objects
      sim.H0$logb.H0 <- weibull(sim.H0$n, a = a.H0, asympt = asympt.H0, lambda = lambda, k = k)
      sim.H0$logb.H1 <- weibull(sim.H0$n, a = a.H1, asympt = asympt.H1, lambda = lambda, k = k)
      
      
      sim.H1$logb.H0 <- weibull(sim.H1$n, a = a.H0, asympt = asympt.H0, lambda = lambda, k = k)
      sim.H1$logb.H1 <- weibull(sim.H1$n, a = a.H1, asympt = asympt.H1, lambda = lambda, k = k)
      
    }
    
    # Get Results under H0
    res.H0 <- sim.H0 %>%
      filter(.[[BF.column]] >= logb.H1 | .[[BF.column]] <= logb.H0) %>%
      group_by(id) %>%
      filter(row_number()==1) %>%
      ungroup()
    
    # Get Results under H1
    res.H1 <- sim.H1 %>%
      filter(.[[BF.column]] >= logb.H1 | .[[BF.column]] <= logb.H0) %>%
      group_by(id) %>%
      filter(row_number()==1) %>%
      ungroup()
    
  }
  
  if(method == "collapsing" | method == "collapsing.sym"){
    med.H1 <- median(res.H1$n)
    mean.H1 <- mean(res.H1$n)
    FN.rate <- nrow(res.H1 %>% filter(.[[BF.column]] <= logb.H0))/nrow(res.H1)
    med.H0 <- median(res.H0$n)
    mean.H0 <- mean(res.H0$n)
    FP.rate <- nrow(res.H0 %>% filter(.[[BF.column]] >= logb.H1))/nrow(res.H0)
  } else {
    med.H0 <- median(res.H0$n)
    mean.H0 <- mean(res.H0$n)
    FP.rate <- nrow(res.H0 %>% filter(.[[BF.column]] >= log.a))/nrow(res.H0)
    med.H1 <- median(res.H1$n)
    mean.H1 <- mean(res.H1$n)
    FN.rate <- nrow(res.H1 %>% filter(.[[BF.column]] <= log.b))/nrow(res.H1)
    
  }
  
  results <- list(res.H0 = res.H0, res.H1 = res.H1, median.H1 = med.H1, median.H0 = med.H0, mean.H1 = mean.H1, mean.H0 = mean.H0, FN.rate = FN.rate, FP.rate = FP.rate)
  
  return(results)
  
  
}