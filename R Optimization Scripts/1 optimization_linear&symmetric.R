# ------------------------------------------------------------------------------
# -------------- OPTIMIZATION WITH LINEAR SYMMETRIC BOUNDARIES -----------------
# ------------------------------------------------------------------------------

#'@import dplyr
#'
#'@load minimize.medianN
#'@source seq.simresults (and choose sim.H0 and sim.H1)

# ---------------- Optimization with the base package --------------------------
sim.H0 <- SequentialBF.obsub.0$sim
sim.H1 <- SequentialBF.obsub.0.8$sim

fn <- function(bounds){
  N <- minimize.meanN(bounds,
                   sim.H0 = sim.H0,
                   sim.H1 = sim.H1,
                   max.FP = 0.05,
                   max.FN = 0.05,
                   BF.column = 5)
  return(N)
}

base.opt.mean.DIST <- optimize(f = fn,
                               interval = c(1,30),
                               maximum = FALSE)

