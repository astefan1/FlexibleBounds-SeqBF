x# ------------------------------------------------------------------------------
# ------------- OPTIMIZATION WITH COLLAPSING BOUNDARIES (WEIBULL) --------------
# ------------------------------------------------------------------------------

#'@import dplyr
#'
#'@load minimize.medianN.r
#'@source seq.simresults (and choose sim.H0 and sim.H1)

# ------------ Differential Evolution Optimization with NMOF -------------------

#' @import NMOF

### Condition 5, ES = 0.5 / 0.2 / 0.8 ###

sim.H0 <- seq.simresults$SequentialBF.obsub.0$sim
sim.H1 <- seq.simresults$SequentialBF.obsub.0.5$sim
sim.H1 <- seq.simresults$SequentialBF.obsub.0.2$sim

OF.DE <- function(params){
  N <- minimize.meanN.r(params,
                        sim.H0 = sim.H0,
                        sim.H1 = sim.H1,
                        max.FP = 0.05,
                        max.FN = 0.2,
                        BF.column = 5)
  return(N)
}

# initP <- matrix(c(log(seq(from = 3, to = 30, length.out = 20)),
#                   log(seq(from = 2, to = 29, length.out = 20)),
#                   seq(from = 1, to = 1050, length.out = 20),
#                   seq(from = 0, to = 100, length.out = 20)),
#                 byrow = TRUE,
#                 nrow = 4)

set.seed(2819475)
initP <- matrix(c(runif(20, min = log(1), max = log(30)),
                  runif(20, min = log(1), max = log(30)),
                  runif(20, min = 0, max = 1050),
                  runif(20, min = 0, max = 100)),
     byrow = TRUE,
     nrow = 4)

algo.DE <- list(CR = 0.8,
                F = 0.5,
                nP = 20,
                nG = 200,
                min = c(0, 0, 0, 0),
                max = c(log(30), log(30), 1050, 100),
                minmaxConstr = TRUE,
                initP = initP,
                printDetail = TRUE,
                printBar = TRUE,
                storeF = TRUE,
                storeSolutions = TRUE)

opt.DE <- DEopt(OF.DE, algo.DE)
