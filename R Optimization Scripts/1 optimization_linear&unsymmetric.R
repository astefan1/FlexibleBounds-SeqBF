# ------------------------------------------------------------------------------
# ---------------- OPTIMIZATION WITH LINEAR BOUNDARIES -------------------------
# ------------------------------------------------------------------------------

#'@import dplyr
#'
#'@load minimize.medianN
#'@source seq.simresults (and choose sim.H0 and sim.H1)

# ---------------- Optimization with the NMKB package --------------------------

# This does not work properly because algorithm gets stuck at local minima
# We need a different optimization algorithm

#'@import dfoptim

fn <- function(bounds){
  N <- minimize.medianN(bounds,
                   sim.H0 = sim.H0,
                   sim.H1 = sim.H1,
                   max.FP = 0.05,
                   max.FN = 0.2,
                   BF.column = 5)
  return(N)
}

par = c(1/20,5)

nmkb(par = par,
    fn = fn,
    lower = c(1/30,1),
    upper = c(1,30),
    control = list(maxfeval = 50,
                   tol = 0.001,
                   regsimp = TRUE,
                   maximize = FALSE,
                   restarts.max = 3,
                   trace = FALSE))

# -------------- Particle Swarm Optimization with NMOF -------------------------

#'@import NMOF

OF.PSO <- function(bounds){
  N <- minimize.medianN(bounds,
                        sim.H0 = sim.H0,
                        sim.H1 = sim.H1,
                        max.FP = 0.05,
                        max.FN = 0.2,
                        BF.column = 5)
  return(N)
}

nP <- 10
initP <- matrix(c(1/3, 1/6, 1/9, 1/12, 1/15, 1/18, 1/21, 1/24, 1/27, 1/30,
                  3, 6, 9, 12, 15, 18, 21, 24, 27, 30),
                ncol = 10,
                nrow = 2,
                byrow = TRUE)

algo.PSO <- list(nP = nP,
             nG = 100,
             c1 = 1,
             c2 = 2,
             iner = 0.9,
             initV = 1,
             minmaxConstr = TRUE,
             min = c(1/30,1),
             max = c(1, 30),
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
