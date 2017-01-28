# canprot/R/diffplot.R
# plot mean differences of ZC and nH2O, or other variables
# 20160715 jmd

diffplot <- function(comptab, vars=c("ZC", "nH2O"), col="black", plot.rect=FALSE, plot.text=TRUE) {
  # convert to data frame if needed
  if(!is.data.frame(comptab)) comptab <- do.call(rbind, comptab)
  # which columns we're using
  stats <- c("diff", "CLES", "p.value")
  iX <- sapply(paste(vars[1], stats, sep="."), grep, colnames(comptab))
  iY <- sapply(paste(vars[2], stats, sep="."), grep, colnames(comptab))
  # get mean difference, common language effect size and p-value
  X_d <- comptab[, iX[1]]
  X_e <- signif(comptab[, iX[2]], 2)
  X_p <- comptab[, iX[3]]
  Y_d <- comptab[, iY[1]]
  Y_e <- signif(comptab[, iY[2]], 2)
  Y_p <- comptab[, iY[3]]
  # set up plot
  x <- vars[1]
  y <- vars[2]
  if(vars[1]=="ZC") x <- cplab$ZC[[1]]
  if(vars[2]=="nH2O") y <- cplab$nH2O[[1]]
  xlab <- substitute("mean difference (" * x * ")", list(x=x))
  ylab <- substitute("mean difference (" * y * ")", list(y=y))
  # initialize plot: add a 0 to make sure we can see the axis
  plot(type="n", c(X_d, 0), c(Y_d, 0), xlab=xlab, ylab=ylab)
  # add a reference rectangle
  if(plot.rect) rect(-0.01, -0.01, 0.01, 0.01, col="grey80", lwd=0)
  # show axis lines
  abline(h=0, lty=3)
  abline(v=0, lty=3)
  # show drop lines: dotted/solid if p-value/effect size meet criteria
  lty.X <- ifelse(abs(X_e - 50) >= 10, 1, ifelse(X_p < 0.05, 2, 0))
  lty.Y <- ifelse(abs(Y_e - 50) >= 10, 1, ifelse(Y_p < 0.05, 2, 0))
  for(i in seq_along(X_d)) {
    lines(rep(X_d[i], 2), c(0, Y_d[i]), lty=lty.X[i])
    lines(c(0, X_d[i]), rep(Y_d[i], 2), lty=lty.Y[i])
  } 
  # point symbols: open circle, filled circle, filled square (0, 1 or 2 vars with p-value < 0.05)
  p_signif <- rowSums(data.frame(X_p < 0.05, Y_p < 0.05))
  pch <- ifelse(p_signif==2, 15, ifelse(p_signif==1, 19, 21))
  # plot points with specified color
  col <- rep(col, length.out=nrow(comptab))
  if(!plot.text) points(X_d, Y_d, pch=pch, col=col, bg="white")
  else {
    # plot bigger points with letters inside
    points(X_d, Y_d, pch=pch, col=col, bg="white", cex=2)
    # use white letters on colored background
    col[pch %in% c(15, 19)] <- "white"
    text(X_d, Y_d, c(letters, LETTERS)[seq_along(X_d)], col=col, cex=0.9)
  }
}
