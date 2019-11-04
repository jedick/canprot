# canprot/R/diffplot.R
# plot mean or median differences of ZC and nH2O, or other variables
# 20160715 jmd
# 20190329 add oldstyle = FALSE (no drop lines; show kernel density)

diffplot <- function(comptab, vars=c("ZC", "nH2O"), col="black", plot.rect=FALSE, pt.text=c(letters, LETTERS),
                     cex.text = 0.9, oldstyle = TRUE, pch = 1, cex = 2, contour = TRUE, col.contour = par("fg")) {
  # convert to data frame if needed
  if(!is.data.frame(comptab)) comptab <- do.call(rbind, comptab)
  # which columns we're using
  stats <- c("diff", "CLES", "p.value")
  iX <- sapply(paste(vars[1], stats, sep="."), grep, colnames(comptab))
  iY <- sapply(paste(vars[2], stats, sep="."), grep, colnames(comptab))
  # get mean/median difference, common language effect size and p-value
  X_d <- comptab[, iX[1]]
  X_e <- signif(comptab[, iX[2]], 2)
  X_p <- comptab[, iX[3]]
  Y_d <- comptab[, iY[1]]
  Y_e <- signif(comptab[, iY[2]], 2)
  Y_p <- comptab[, iY[3]]
  # set up plot
  Dx <- paste0("D", vars[1])
  Dy <- paste0("D", vars[2])
  if(oldstyle) {
    x <- cplabbar[[Dx]][[1]]
    y <- cplabbar[[Dy]][[1]]
  } else {
    x <- cplab[[Dx]][[1]]
    y <- cplab[[Dy]][[1]]
  }
  # use colnames to figure out whether the difference is of the mean or median
  if(any(grepl("mean", colnames(comptab)))) mfun <- "mean"
  if(any(grepl("median", colnames(comptab)))) mfun <- "median"
  xlab <- substitute(mfun * " difference (" * x * ")", list(mfun=mfun, x=x))
  ylab <- substitute(mfun * " difference (" * y * ")", list(mfun=mfun, y=y))
  # initialize plot: add a 0 to make sure we can see the axis
  plot(type="n", c(X_d, 0), c(Y_d, 0), xlab=xlab, ylab=ylab)
  # contour 2-D kernel density estimate 20190329
  # https://stats.stackexchange.com/questions/31726/scatterplot-with-contour-heat-overlay
  if(!oldstyle & any(contour)) {
    # include points specified by 'contour' 20191102
    z <- kde2d(X_d[contour], Y_d[contour], n = 50)
    contour(z, drawlabels = FALSE, nlevels = 3, lty = 2, add = TRUE, col = col.contour)
  }
  # add a reference rectangle
  if(plot.rect) rect(-0.01, -0.01, 0.01, 0.01, col="grey80", lwd=0)
  # show axis lines
  abline(h=0, lty=3)
  abline(v=0, lty=3)
  if(oldstyle) {
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
  }
  # plot points with specified color (and point style, only for oldstyle = FALSE)
  col <- rep(col, length.out=nrow(comptab))
  pch <- rep(pch, length.out=nrow(comptab))
  points(X_d, Y_d, pch=pch, col=col, bg="white", cex=cex)
  if(!identical(pt.text, NA) | !identical(pt.text, FALSE)) {
    # add letters; for dark background, use white letters
    col[pch %in% c(15, 19)] <- "white"
    text(X_d, Y_d, pt.text[seq_along(X_d)], col=col, cex=cex.text)
  }
}
