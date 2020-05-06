# canprot/qdist.R
# 20200505 jmd

# quantile distribution for up- and down-regulated proteins
# first version 20200428
# use ecdf() to calculate knots 20200506
qdist <- function(pd = pdat_colorectal("JKMF10"), vars = c("ZC", "nH2O"), show.steps = FALSE) {
  # initialize plot
  if(length(vars)==2) par(mfrow = c(2, 1))
  par(yaxs = "i")
  for(var in vars) {
    if(var=="ZC") {
      up <- pd$pcomp$ZC[pd$up2]
      dn <- pd$pcomp$ZC[!pd$up2]
      xlab <- cplab$ZC
    }
    if(var=="nH2O") {
      up <- pd$pcomp$residue.basis[, "H2O"][pd$up2]
      dn <- pd$pcomp$residue.basis[, "H2O"][!pd$up2]
      xlab <- cplab$nH2O
    }
    # start the plot
    plot(range(up, dn), c(0, 1), type = "n", xlab = xlab, ylab = "Quantile point")
    # a function to plot the points and lines
    pfun <- function(x, ...) {
      Fn <- ecdf(x)
      # plot the values (knots) and verticals
      if(show.steps) plot(Fn, add = TRUE, cex = 0.25, col.01line = NA, verticals = TRUE, col.hor = "gray70", ...)
      # add lines that bisect the verticals (to intersect the quantile points)
      x <- knots(Fn)
      y <- sort(Fn(x)) - 0.5 / length(x)
      lines(x, y, ...)
    }
    pfun(dn)
    pfun(up, col = 2, lty = 2)
    # draw a line at the 0.5 quantile
    lines(c(median(dn), median(up)), c(0.5, 0.5), lwd = 2, col = "slategray3")
  }
}
