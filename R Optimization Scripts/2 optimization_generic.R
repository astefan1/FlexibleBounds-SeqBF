# ------------------------------------------------------------------------------
# ------- VARIABLE OPTIMIZATION FUNCTION FOR BOUNDARY OPTIMIZATION -------------
# ------------------------------------------------------------------------------

#'@source 1 minimize.meanN_linear&symmetric
#'@source 1 minimize.medianN_linear&symmetric.R
#'@source 1 minimize.meanN_linear&unsymmetric
#'@source 1 minimize.medianN_linear&unsymmetric
#'@source 2 minimize.meanN_collapsing restrained
#'@source 2 minimize.medianN_collapsing restrained
#'@source 2 collapsing bounds function weibull
#'
#'@import NMOF
#'@import dplyr
#'

bounds.opt <- function(sim.H0, sim.H1, max.FP, max.FN, initP=NULL, type = "linear", criterion = "mean", BF.column=5, nP=20, nG=200){
  
  # Define arguments
  match.arg(initP, c("set1", "set2", NULL), several.ok = FALSE)
  match.arg(type, c("linear", "collapsing.sym", "linear.sym"), several.ok = FALSE)
  match.arg(criterion, c("mean", "median"), several.ok = FALSE)
  
  # Define objective function
  
  shape <- ifelse(type == "linear", "linear", ifelse(type == "collapsing.sym", "r", "s"))
  minimize <- get(paste0("minimize.", criterion, "N.", shape))
  
  OF.DE <- function(params){
    N <- minimize(params,
                    sim.H0 = sim.H0,
                    sim.H1 = sim.H1,
                    max.FP = max.FP,
                    max.FN = max.FN,
                  BF.column = BF.column)
    return(N)
  }
  
  # Define initial parameters for DE algorithm
  if(type == "linear"){
    if(initP == "set1"){
      initP <- matrix(c(1/seq(3,30, length.out = nP),
                        seq(3,30, length.out = nP)),
                      ncol = 10,
                      nrow = 2,
                      byrow = TRUE)
    } else {
      set.seed(2819475)
      initP <- matrix(c(runif(nP, min = 1/30, max = 1),
                        runif(nP, min = 1, max = 30)),
                     ncol = 10,
                     nrow = 2,
                     byrow = TRUE)
    }
  } else if (type == "collapsing.sym"){
    if(initP == "set1"){
      initP <- matrix(c(log(seq(from = 3, to = 30, length.out = nP)),
                        log(seq(from = 2, to = 29, length.out = nP)),
                        seq(from = 1, to = 1050, length.out = nP),
                        seq(from = 0, to = 100, length.out = nP)),
                      byrow = TRUE,
                      nrow = 4)
    } else{
      set.seed(2819475)
      initP <- matrix(c(runif(nP, min = log(1), max = log(30)),
                        runif(nP, min = log(1), max = log(30)),
                        runif(nP, min = 0, max = 1050),
                        runif(nP, min = 0, max = 100)),
                      byrow = TRUE,
                      nrow = 4)
    }
  } 
  
  # Optimization algorithm specification
  if(type == "linear"){
    minima <- c(1/30, 1) 
    maxima <- c(1, 30)
  } else if (type == "collapsing.sym"){
    minima <- c(0, 0, 0, 0)
    maxima <- c(log(30), log(30), 1050, 100)
  }
  
  if(type == "linear.sym"){
    opt.DE <- optimize(f = OF.DE,
                                   interval = c(1,30),
                                   maximum = FALSE)
  } else {
    algo.DE <- list(CR = 0.8,
                    F = 0.5,
                    nP = nP,
                    nG = nG,
                    min = minima,
                    max = maxima,
                    minmaxConstr = TRUE,
                    initP = initP,
                    printDetail = TRUE,
                    printBar = TRUE,
                    storeF = TRUE,
                    storeSolutions = TRUE)
    
    # Optimization Execution
    opt.DE <- NMOF::DEopt(OF.DE, algo.DE)
  }
 
 return(opt.DE)
  
  }
