weibull <- function(t, a, asympt, lambda, k){
  y <- a - (1 - exp(1)**(-(t/lambda)**k))*(a-asympt)
  return(y)
}

# curve(weibull(x, 12, 3, 0.5, 3), xlim = c(0, 100))


