setwd("~/Documents")

source("2 optimization_generic.R")
source("2 minimize.medianN_collapsing restrained.R")
source("2 minimize.meanN_collapsing restrained.R")
source("2 collapsing bounds function weibull.R")
source("1 minimize.medianN_linear&unsymmetric.R")
source("1 minimize.meanN_linear&unsymmetric.R")
source("1 minimize.medianN_linear&symmetric.R")
source("1 minimize.meanN_linear&symmetric.R")

load("SequentialBF.obsub.0.RData")
load("SequentialBF.obsub.0.2.RData")
load("SequentialBF.obsub.0.5.RData")
load("SequentialBF.obsub.0.8.RData")
load("SequentialBF.obsub.tdist.RData")
sim.H0 <- SequentialBF.obsub.0$sim
sim.H1.02 <- SequentialBF.obsub.0.2$sim
sim.H1.05 <- SequentialBF.obsub.0.5$sim
sim.H1.08 <- SequentialBF.obsub.0.8$sim
sim.H1.DIST <- SequentialBF.obsub.tdist$sim

# Constant, symmetric: ES = 0.2

base.opt.0.2 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.02, max.FP = 0.05, max.FN = 0.2, type = "linear.sym", criterion = "median", BF.column = 5)
base.opt.0.2.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.02, max.FP = 0.05, max.FN = 0.05, type = "linear.sym", criterion = "median", BF.column = 5)
base.opt.mean.0.2 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.02, max.FP = 0.05, max.FN = 0.2, type = "linear.sym", criterion = "mean", BF.column = 5)
base.opt.mean.0.2.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.02, max.FP = 0.05, max.FN = 0.05, type = "linear.sym", criterion = "mean", BF.column = 5)

# Constant, symmetric: ES = 0.5

base.opt.0.5 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.05, max.FP = 0.05, max.FN = 0.2, type = "linear.sym", criterion = "median", BF.column = 5)
base.opt.0.5.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.05, max.FP = 0.05, max.FN = 0.05, type = "linear.sym", criterion = "median", BF.column = 5)
base.opt.mean.0.5 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.05, max.FP = 0.05, max.FN = 0.2, type = "linear.sym", criterion = "mean", BF.column = 5)
base.opt.mean.0.5.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.05, max.FP = 0.05, max.FN = 0.05, type = "linear.sym", criterion = "mean", BF.column = 5)

# Constant, symmetric: ES = 0.8
base.opt.0.8 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.08, max.FP = 0.05, max.FN = 0.2, type = "linear.sym", criterion = "median", BF.column = 5)
base.opt.0.8.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.08, max.FP = 0.05, max.FN = 0.05, type = "linear.sym", criterion = "median", BF.column = 5)
base.opt.mean.0.8 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.08, max.FP = 0.05, max.FN = 0.2, type = "linear.sym", criterion = "mean", BF.column = 5)
base.opt.mean.0.8.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.08, max.FP = 0.05, max.FN = 0.05, type = "linear.sym", criterion = "mean", BF.column = 5)

# Constant, symmetric: ES = DIST
base.opt.DIST <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.DIST, max.FP = 0.05, max.FN = 0.2, type = "linear.sym", criterion = "median", BF.column = 5)
base.opt.DIST.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.DIST, max.FP = 0.05, max.FN = 0.05, type = "linear.sym", criterion = "median", BF.column = 5)
base.opt.mean.DIST <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.DIST, max.FP = 0.05, max.FN = 0.2, type = "linear.sym", criterion = "mean", BF.column = 5)
base.opt.mean.DIST.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.DIST, max.FP = 0.05, max.FN = 0.05, type = "linear.sym", criterion = "mean", BF.column = 5)

# Constant, non-symmetric: ES = 0.2

DE.opt.linear.1.02 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.02, initP = "set1", max.FP = 0.05, max.FN = 0.2, type = "linear", criterion = "median", BF.column = 5, nP = 10, nG = 100)
DE.opt.linear.2.02 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.02, initP = "set2", max.FP = 0.05, max.FN = 0.2, type = "linear", criterion = "median", BF.column = 5, nP = 10, nG = 100)
DE.opt.linear.1.02.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.02, initP = "set1", max.FP = 0.05, max.FN = 0.05, type = "linear", criterion = "median", BF.column = 5, nP = 10, nG = 100)
DE.opt.linear.2.02.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.02, initP = "set2", max.FP = 0.05, max.FN = 0.05, type = "linear", criterion = "median", BF.column = 5, nP = 10, nG = 100)
opt.DE.means.linear.0.2 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.02, initP = "set1", max.FP = 0.05, max.FN = 0.2, type = "linear", criterion = "mean", BF.column = 5, nP = 10, nG = 100)
opt.DE.means.linear_2_.0.2 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.02, initP = "set2", max.FP = 0.05, max.FN = 0.2, type = "linear", criterion = "mean", BF.column = 5, nP = 10, nG = 100)
opt.DE.means.linear.0.2.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.02, initP = "set1", max.FP = 0.05, max.FN = 0.05, type = "linear", criterion = "mean", BF.column = 5, nP = 10, nG = 100)
opt.DE.means.linear_2_.0.2.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.02, initP = "set2", max.FP = 0.05, max.FN = 0.05, type = "linear", criterion = "mean", BF.column = 5, nP = 10, nG = 100)

# Constant, non-symmetric: ES = 0.5

DE.opt.linear.1.05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.05, initP = "set1", max.FP = 0.05, max.FN = 0.2, type = "linear", criterion = "median", BF.column = 5, nP = 10, nG = 100)
DE.opt.linear.2.05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.05, initP = "set2", max.FP = 0.05, max.FN = 0.2, type = "linear", criterion = "median", BF.column = 5, nP = 10, nG = 100)
DE.opt.linear.1.05.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.05, initP = "set1", max.FP = 0.05, max.FN = 0.05, type = "linear", criterion = "median", BF.column = 5, nP = 10, nG = 100)
DE.opt.linear.2.05.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.05, initP = "set2", max.FP = 0.05, max.FN = 0.05, type = "linear", criterion = "median", BF.column = 5, nP = 10, nG = 100)
opt.DE.means.linear.0.5 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.05, initP = "set1", max.FP = 0.05, max.FN = 0.2, type = "linear", criterion = "mean", BF.column = 5, nP = 10, nG = 100)
opt.DE.means.linear_2_.0.5 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.05, initP = "set2", max.FP = 0.05, max.FN = 0.2, type = "linear", criterion = "mean", BF.column = 5, nP = 10, nG = 100)
opt.DE.means.linear.0.5.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.05, initP = "set1", max.FP = 0.05, max.FN = 0.05, type = "linear", criterion = "mean", BF.column = 5, nP = 10, nG = 100)
opt.DE.means.linear_2_.0.5.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.05, initP = "set2", max.FP = 0.05, max.FN = 0.05, type = "linear", criterion = "mean", BF.column = 5, nP = 10, nG = 100)

# Constant, non-symmetric: ES = 0.8

DE.opt.linear.1.08 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.08, initP = "set1", max.FP = 0.05, max.FN = 0.2, type = "linear", criterion = "median", BF.column = 5, nP = 10, nG = 100)
DE.opt.linear.2.08 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.08, initP = "set2", max.FP = 0.05, max.FN = 0.2, type = "linear", criterion = "median", BF.column = 5, nP = 10, nG = 100)
DE.opt.linear.1.08.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.08, initP = "set1", max.FP = 0.05, max.FN = 0.05, type = "linear", criterion = "median", BF.column = 5, nP = 10, nG = 100)
DE.opt.linear.2.08.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.08, initP = "set2", max.FP = 0.05, max.FN = 0.05, type = "linear", criterion = "median", BF.column = 5, nP = 10, nG = 100)
opt.DE.means.linear.0.8 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.08, initP = "set1", max.FP = 0.05, max.FN = 0.2, type = "linear", criterion = "mean", BF.column = 5, nP = 10, nG = 100)
opt.DE.means.linear_2_.0.8 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.08, initP = "set2", max.FP = 0.05, max.FN = 0.2, type = "linear", criterion = "mean", BF.column = 5, nP = 10, nG = 100)
opt.DE.means.linear.0.8.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.08, initP = "set1", max.FP = 0.05, max.FN = 0.05, type = "linear", criterion = "mean", BF.column = 5, nP = 10, nG = 100)
opt.DE.means.linear_2_.0.8.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.08, initP = "set2", max.FP = 0.05, max.FN = 0.05, type = "linear", criterion = "mean", BF.column = 5, nP = 10, nG = 100)

# Constant, non-symmetric: ES = t-dist

DE.opt.linear.1.DIST <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.DIST, initP = "set1", max.FP = 0.05, max.FN = 0.2, type = "linear", criterion = "median", BF.column = 5, nP = 10, nG = 100)
DE.opt.linear.2.DIST <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.DIST, initP = "set2", max.FP = 0.05, max.FN = 0.2, type = "linear", criterion = "median", BF.column = 5, nP = 10, nG = 100)
DE.opt.linear.1.DIST.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.DIST, initP = "set1", max.FP = 0.05, max.FN = 0.05, type = "linear", criterion = "median", BF.column = 5, nP = 10, nG = 100)
DE.opt.linear.2.DIST.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.DIST, initP = "set2", max.FP = 0.05, max.FN = 0.05, type = "linear", criterion = "median", BF.column = 5, nP = 10, nG = 100)
opt.DE.means.linear.DIST <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.DIST, initP = "set1", max.FP = 0.05, max.FN = 0.2, type = "linear", criterion = "mean", BF.column = 5, nP = 10, nG = 100)
opt.DE.means.linear_2_.DIST <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.DIST, initP = "set2", max.FP = 0.05, max.FN = 0.2, type = "linear", criterion = "mean", BF.column = 5, nP = 10, nG = 100)
opt.DE.means.linear.DIST.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.DIST, initP = "set1", max.FP = 0.05, max.FN = 0.05, type = "linear", criterion = "mean", BF.column = 5, nP = 10, nG = 100)
opt.DE.means.linear_2_.DIST.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.DIST, initP = "set2", max.FP = 0.05, max.FN = 0.05, type = "linear", criterion = "mean", BF.column = 5, nP = 10, nG = 100)

# Collapsing, symmetric: ES = 0.2
DE.opt.r.1.0.2 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.02, initP = "set1", max.FP = 0.05, max.FN = 0.2, type = "collapsing.sym", criterion = "median", BF.column = 5, nP = 20, nG = 200)
DE.opt.r.2.0.2 <- bounds.opt (sim.H0 = sim.H0, sim.H1 = sim.H1.02, initP = "set2", max.FP = 0.05, max.FN = 0.2, type = "collapsing.sym", criterion = "median", BF.column = 5, nP = 20, nG = 200)
DE.opt.r.1.0.2.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.02, initP = "set1", max.FP = 0.05, max.FN = 0.05, type = "collapsing.sym", criterion = "median", BF.column = 5, nP = 20, nG = 200)
DE.opt.r.2.0.2.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.02, initP = "set2", max.FP = 0.05, max.FN = 0.05, type = "collapsing.sym", criterion = "median", BF.column = 5, nP = 20, nG = 200)
# opt.DE.means.r.0.2 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.02, initP = "set1", max.FP = 0.05, max.FN = 0.2, type = "collapsing.sym", criterion = "mean", BF.column = 5, nP = 20, nG = 200)
# opt.DE.means.r.0.2_2 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.02, initP = "set2", max.FP = 0.05, max.FN = 0.2, type = "collapsing.sym", criterion = "mean", BF.column = 5, nP = 20, nG = 200)
opt.DE.means.r.0.2.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.02, initP = "set1", max.FP = 0.05, max.FN = 0.05, type = "collapsing.sym", criterion = "mean", BF.column = 5, nP = 20, nG = 200)
opt.DE.means.r.0.2_2.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.02, initP = "set2", max.FP = 0.05, max.FN = 0.05, type = "collapsing.sym", criterion = "mean", BF.column = 5, nP = 20, nG = 200)

# Collapsing, symmetric: ES = 0.5
# DE.opt.r.1.0.5 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.05, initP = "set1", max.FP = 0.05, max.FN = 0.2, type = "collapsing.sym", criterion = "median", BF.column = 5, nP = 20, nG = 200)
# DE.opt.r.2.0.5 <- bounds.opt (sim.H0 = sim.H0, sim.H1 = sim.H1.05, initP = "set2", max.FP = 0.05, max.FN = 0.2, type = "collapsing.sym", criterion = "median", BF.column = 5, nP = 20, nG = 200)
DE.opt.r.1.0.5.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.05, initP = "set1", max.FP = 0.05, max.FN = 0.05, type = "collapsing.sym", criterion = "median", BF.column = 5, nP = 20, nG = 200)
DE.opt.r.2.0.5.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.05, initP = "set2", max.FP = 0.05, max.FN = 0.05, type = "collapsing.sym", criterion = "median", BF.column = 5, nP = 20, nG = 200)
# opt.DE.means.r.0.5 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.05, initP = "set1", max.FP = 0.05, max.FN = 0.2, type = "collapsing.sym", criterion = "mean", BF.column = 5, nP = 20, nG = 200)
# opt.DE.means.r.0.5_2 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.05, initP = "set2", max.FP = 0.05, max.FN = 0.2, type = "collapsing.sym", criterion = "mean", BF.column = 5, nP = 20, nG = 200)
# opt.DE.means.r.0.5.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.05, initP = "set1", max.FP = 0.05, max.FN = 0.05, type = "collapsing.sym", criterion = "mean", BF.column = 5, nP = 20, nG = 200)
# opt.DE.means.r.0.5_2.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.05, initP = "set2", max.FP = 0.05, max.FN = 0.05, type = "collapsing.sym", criterion = "mean", BF.column = 5, nP = 20, nG = 200)

# Collapsing, symmetric: ES = 0.8
DE.opt.r.1.0.8 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.08, initP = "set1", max.FP = 0.05, max.FN = 0.2, type = "collapsing.sym", criterion = "median", BF.column = 5, nP = 20, nG = 200)
DE.opt.r.2.0.8 <- bounds.opt (sim.H0 = sim.H0, sim.H1 = sim.H1.08, initP = "set2", max.FP = 0.05, max.FN = 0.2, type = "collapsing.sym", criterion = "median", BF.column = 5, nP = 20, nG = 200)
save(DE.opt.r.2.0.8, file = "DE.opt.r.2.0.8.RData")
DE.opt.r.1.0.8.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.08, initP = "set1", max.FP = 0.05, max.FN = 0.05, type = "collapsing.sym", criterion = "median", BF.column = 5, nP = 20, nG = 200)
save(DE.opt.r.1.0.8.ME05, file = "DE.opt.r.1.0.8.ME05.RData")
DE.opt.r.2.0.8.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.08, initP = "set2", max.FP = 0.05, max.FN = 0.05, type = "collapsing.sym", criterion = "median", BF.column = 5, nP = 20, nG = 200)
save(DE.opt.r.2.0.8.ME05, file = "DE.opt.r.2.0.8.ME05.RData")
opt.DE.means.r.0.8 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.08, initP = "set1", max.FP = 0.05, max.FN = 0.2, type = "collapsing.sym", criterion = "mean", BF.column = 5, nP = 20, nG = 200)
save(opt.DE.means.r.0.8, file = "opt.DE.means.r.0.8.RData")
opt.DE.means.r.0.8_2 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.08, initP = "set2", max.FP = 0.05, max.FN = 0.2, type = "collapsing.sym", criterion = "mean", BF.column = 5, nP = 20, nG = 200)
save(opt.DE.means.r.0.8_2, file = "opt.DE.means.r.0.8_2.RData")
opt.DE.means.r.0.8.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.08, initP = "set1", max.FP = 0.05, max.FN = 0.05, type = "collapsing.sym", criterion = "mean", BF.column = 5, nP = 20, nG = 200)
save(opt.DE.means.r.0.8.ME05, file = "opt.DE.means.r.0.8.ME05.RData")
opt.DE.means.r.0.8_2.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.08, initP = "set2", max.FP = 0.05, max.FN = 0.05, type = "collapsing.sym", criterion = "mean", BF.column = 5, nP = 20, nG = 200)
save(opt.DE.means.r.0.8_2.ME05, file = "opt.DE.means.r.0.8_2.ME05.RData")

# Collapsing, symmetric: ES = t-dist
DE.opt.r.1.DIST <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.DIST, initP = "set1", max.FP = 0.05, max.FN = 0.2, type = "collapsing.sym", criterion = "median", BF.column = 5, nP = 20, nG = 200)
save(DE.opt.r.1.DIST, file = "DE.opt.r.1.DIST.RData")
DE.opt.r.2.DIST <- bounds.opt (sim.H0 = sim.H0, sim.H1 = sim.H1.DIST, initP = "set2", max.FP = 0.05, max.FN = 0.2, type = "collapsing.sym", criterion = "median", BF.column = 5, nP = 20, nG = 200)
save(DE.opt.r.2.DIST, file = "DE.opt.r.2.DIST.RData")
DE.opt.r.1.DIST.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.DIST, initP = "set1", max.FP = 0.05, max.FN = 0.05, type = "collapsing.sym", criterion = "median", BF.column = 5, nP = 20, nG = 200)
save(DE.opt.r.1.DIST.ME05, file = "DE.opt.r.1.DIST.ME05.RData")
DE.opt.r.2.DIST.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.DIST, initP = "set2", max.FP = 0.05, max.FN = 0.05, type = "collapsing.sym", criterion = "median", BF.column = 5, nP = 20, nG = 200)
save(DE.opt.r.2.DIST.ME05, file = "DE.opt.r.2.DIST.ME05.RData")
opt.DE.means.r.DIST <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.DIST, initP = "set1", max.FP = 0.05, max.FN = 0.2, type = "collapsing.sym", criterion = "mean", BF.column = 5, nP = 20, nG = 200)
save(opt.DE.means.r.DIST, file = "opt.DE.means.r.DIST.RData")
opt.DE.means.r.DIST_2 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.DIST, initP = "set2", max.FP = 0.05, max.FN = 0.2, type = "collapsing.sym", criterion = "mean", BF.column = 5, nP = 20, nG = 200)
save(opt.DE.means.r.DIST_2, file = "opt.DE.means.r.DIST_2.RData")
opt.DE.means.r.DIST.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.DIST, initP = "set1", max.FP = 0.05, max.FN = 0.05, type = "collapsing.sym", criterion = "mean", BF.column = 5, nP = 20, nG = 200)
save(opt.DE.means.r.DIST.ME05, file = "opt.DE.means.r.DIST.ME05.RData")
opt.DE.means.r.DIST_2.ME05 <- bounds.opt(sim.H0 = sim.H0, sim.H1 = sim.H1.DIST, initP = "set2", max.FP = 0.05, max.FN = 0.05, type = "collapsing.sym", criterion = "mean", BF.column = 5, nP = 20, nG = 200)
save(opt.DE.means.r.DIST_2.ME05, file = "opt.DE.means.r.DIST_2.ME05.RData")