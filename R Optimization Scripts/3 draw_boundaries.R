draw.boundaries <- function(res.opt, type = c("linear", "collapsing.sym", "symmetric"), title = "", xlim = c(0, 1000), at = seq(0, 1000, by = 200), ...){
  
  if(type == "collapsing.sym"){
    curve(weibull(x, res.opt$xbest[1], res.opt$xbest[2], res.opt$xbest[3], res.opt$xbest[4]), xlim = xlim, ylim = c(log(1/30), log(30)), xaxt = "n", yaxt = "n", xlab = "N", ylab = "BF", bty = "n", ...)
    curve(weibull(x, -res.opt$xbest[1], -res.opt$xbest[2], res.opt$xbest[3], res.opt$xbest[4]), xlim = xlim, ylim = c(log(1/30), log(30)), xaxt = "n", yaxt = "n", xlab = "N", ylab = "BF", bty = "n", add = TRUE, ...)
    axis(side = 1, at = at)
    axis(side = 2, log(c(1/30, 1/10, 1/3, 1, 3, 10, 30)), labels = c("1/30", "1/10", "1/3", "1", "3", "10", "30"), las = 2)
    abline(h = 0, lty = "dotted")
    mtext(title, side = 3)
    
  } else if (type == "linear") {
    plot(NA, xlim = xlim, ylim = c(log(1/30), log(30)), xaxt = "n", yaxt = "n", xlab = "N", ylab = "BF", bty = "n")
    abline(h = 0, lty = "dotted")
    abline(h = log(res.opt$xbest[1]), lty = "solid", ...)
    abline(h = log(res.opt$xbest[2]), lty = "solid", ...)
    axis(side = 1, at = at)
    axis(side = 2, log(c(1/30, 1/10, 1/3, 1, 3, 10, 30)), labels = c("1/30", "1/10", "1/3", "1", "3", "10", "30"), las = 2)
    mtext(title, side = 3)
  } else {
    res.opt$minimum <- c(res.opt$minimum, 1/res.opt$minimum)
    plot(NA, xlim = xlim, ylim = c(log(1/30), log(30)), xaxt = "n", yaxt = "n", xlab = "N", ylab = "BF", bty = "n")
    abline(h = 0, lty = "dotted")
    abline(h = log(res.opt$minimum[1]), lty = "solid", ...)
    abline(h = log(res.opt$minimum[2]), lty = "solid", ...)
    axis(side = 1, at = at)
    axis(side = 2, log(c(1/30, 1/10, 1/3, 1, 3, 10, 30)), labels = c("1/30", "1/10", "1/3", "1", "3", "10", "30"), las = 2)
    mtext(title, side = 3)
  }
}