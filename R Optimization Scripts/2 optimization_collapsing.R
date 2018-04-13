# ------------------------------------------------------------------------------
# ------------- OPTIMIZATION WITH COLLAPSING BOUNDARIES (NORMAL) ---------------
# ------------------------------------------------------------------------------

#'@import dplyr
#'
#'@load minimize.medianN
#'@source seq.simresults (and choose sim.H0 and sim.H1)

# ---------------- Optimization with the NMKB package --------------------------

# This does not work properly because algorithm gets stuck at local minima
# We need a different optimization algorithm

#'@import dfoptim

fn <- function(params){
  sigma.H0 <- params[1]
  sigma.H1 <- params[2]
  k.H0 <- params[3]
  k.H1 <- params[4]
  
  N <- minimize.medianN(sigma.H0 = sigma.H0,
                        sigma.H1 = sigma.H1,
                        k.H0 = k.H0,
                        k.H1 = k.H1,
                        sim.H0 = sim.H0,
                        sim.H1 = sim.H1,
                        max.FP = 0.05,
                        max.FN = 0.2,
                        BF.column = 5)
  return(N)
}

par = c(400,400, 1/10, 20)

try(nmkb(par = par,
    fn = fn,
    lower = c(0.0001, 0.0001, 0, 0),
    upper = c(Inf, Inf, 29/30, 29),
    control = list(maxfeval = 2,
                   tol = 0.001,
                   regsimp = TRUE,
                   maximize = FALSE,
                   restarts.max = 3)))

# -------------- Particle Swarm Optimization with NMOF -------------------------

#'@import NMOF

OF.PSO <- function(params){
  
  sigma.H0 <- params[1]
  sigma.H1 <- params[2]
  k.H0 <- params[3]
  k.H1 <- params[4]
  
  N <- minimize.medianN(sigma.H0 = sigma.H0,
                        sigma.H1 = sigma.H1,
                        k.H0 = k.H0,
                        k.H1 = k.H1,
                        sim.H0 = sim.H0,
                        sim.H1 = sim.H1,
                        max.FP = 0.05,
                        max.FN = 0.2,
                        BF.column = 5)
  return(N)
}

nP <- 10

algo.PSO <- list(nP = nP,
             nG = 100,
             c1 = 1,
             c2 = 2,
             iner = 0.9,
             initV = 1,
             minmaxConstr = TRUE,
             min = c(1, 1, 0, 0),
             max = c(Inf, Inf, 29/30, 30),
             storeSolutions = TRUE,
             storeF = TRUE)

opt <- PSopt(OF.PSO, algo = algo.PSO)

# ------------ Differential Evolution Optimization with NMOF -------------------

#' @import NMOF

OF.DE <- function(bounds){
  N <- minimize.medianN(bounds,
                        sim.H0 = sim.H0,
                        sim.H1 = sim.H1,
                        max.FP = 0.05,
                        max.FN = 0.2,
                        BF.column = 5)
  return(N)
}

initP <- matrix(c(1/3, 1/6, 1/9, 1/12, 1/15, 1/18, 1/21, 1/24, 1/27, 1/30,
                  3, 6, 9, 12, 15, 18, 21, 24, 27, 30),
                ncol = 10,
                nrow = 2,
                byrow = TRUE)

algo.DE <- list(CR = 0.8,
                F = 0.5,
                nP = 10,
                nG = 100,
                min = c(1/30, 1),
                max = c(1, 30),
                minmaxConstr = TRUE,
                initP = initP,
                printDetail = TRUE,
                printBar = TRUE,
                storeF = TRUE,
                storeSolutions = TRUE)

opt.DE <- DEopt(OF.DE, algo.DE)
