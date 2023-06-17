# canprot/qdist.R
# 20200505 jmd

# Quantile distribution for up- and down-regulated proteins
# First version 20200428
# Use ecdf() to calculate knots 20200506
qdist <- function(pdat = pdat_colorectal("JKMF10"), vars = c("Zc", "nH2O"), show.steps = FALSE) {
  # Initialize plot
  if(length(vars)==2) par(mfrow = c(2, 1))
  par(yaxs = "i")
  for(var in vars) {
    if(var=="Zc") {
      X <- Zc(pdat$pcomp$aa)
      xlab <- cplab$Zc
    }
    if(var=="nH2O") {
      X <- nH2O(pdat$pcomp$aa)
      xlab <- cplab$nH2O
    }
    up <- X[pdat$up2]
    dn <- X[!pdat$up2]
    # Start the plot
    plot(range(up, dn), c(0, 1), type = "n", xlab = xlab, ylab = "Quantile point")
    # A function to plot the points and lines
    pfun <- function(x, ...) {
      Fn <- ecdf(x)
      # Plot the values (knots) and verticals
      if(show.steps) plot(Fn, add = TRUE, cex = 0.25, col.01line = NA, verticals = TRUE, col.hor = "gray70", ...)
      # Add lines that bisect the verticals (to intersect the quantile points)
      x <- knots(Fn)
      y <- sort(Fn(x)) - 0.5 / length(x)
      lines(x, y, ...)
    }
    pfun(dn)
    pfun(up, col = 2, lty = 2)
    # Draw a line at the 0.5 quantile
    lines(c(median(dn), median(up)), c(0.5, 0.5), lwd = 2, col = "slategray3")
  }
}
