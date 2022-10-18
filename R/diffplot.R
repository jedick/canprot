# canprot/R/diffplot.R
# Plot mean or median differences of ZC and nH2O, or other variables
# 20160715 jmd
# 20190329 add oldstyle = FALSE (no drop lines; show kernel density)

diffplot <- function(comptab, vars=c("ZC", "nH2O"), col="black", plot.rect=FALSE, pt.text=c(letters, LETTERS),
                     cex.text = 0.85, oldstyle = FALSE, pch = 1, cex = 2.1, contour = TRUE, col.contour = par("fg"),
                     probs = 0.5, add = FALSE, labtext = NULL, ...) {

  # Convert to data frame if needed
  if(!is.data.frame(comptab)) comptab <- do.call(rbind, comptab)
  # Which columns we're using
  stats <- c("diff", "CLES", "p.value")
  iX <- unlist(sapply(paste(vars[1], stats, sep=".*"), grep, colnames(comptab)))
  iY <- unlist(sapply(paste(vars[2], stats, sep=".*"), grep, colnames(comptab)))
  # Get mean/median difference, common language effect size and p-value
  X_d <- comptab[, iX[1]]
  Y_d <- comptab[, iY[1]]

  # Figure out axis labels
  # Only use part before underscore 20191207
  Dx <- paste0("D", strsplit(vars[1], "_")[[1]][1])
  Dy <- paste0("D", strsplit(vars[2], "_")[[1]][1])
  if(oldstyle) {
    # "oldstyle" labels including overbar
    cplabbar <- cplab
    cplabbar$nH2O <- expression(bar(italic(n))[H[2]*O])
    cplabbar$DnH2O <- expression(Delta*bar(italic(n))[H[2]*O])
    xvar <- cplabbar[[Dx]][[1]]
    yvar <- cplabbar[[Dy]][[1]]
    # For oldstyle plots, also get common language effect size and p-value
    X_e <- signif(comptab[, iX[2]], 2)
    X_p <- comptab[, iX[3]]
    Y_e <- signif(comptab[, iY[2]], 2)
    Y_p <- comptab[, iY[3]]
  } else {
    xvar <- cplab[[Dx]][[1]]
    yvar <- cplab[[Dy]][[1]]
  }
  # Use colnames to figure out whether the difference is of the mean or median
  if(is.null(labtext)) {
    # Treat the x- and y-variables separately in case one is median and one is mean (possible with phylostrata) 20191127
    xfun <- gsub("1", "", strsplit(grep(vars[1], colnames(comptab), value = TRUE)[1], "\\.")[[1]][2])
    yfun <- gsub("1", "", strsplit(grep(vars[2], colnames(comptab), value = TRUE)[1], "\\.")[[1]][2])
    xparen <- paste0("(", xfun, " difference)")
    yparen <- paste0("(", yfun, " difference)")
  } else {
    labtext <- rep(labtext, length.out = 2)
    if(is.list(labtext)) lt1 <- labtext[[1]] else lt1 <- labtext[1]
    if(is.list(labtext)) lt2 <- labtext[[2]] else lt2 <- labtext[2]
    xparen <- bquote("("*.(lt1)*")")
    yparen <- bquote("("*.(lt2)*")")
  }
  if(identical(labtext[1], NA)) xlab <- xvar else xlab <- substitute(x ~ xparen, list(xparen=xparen, x=xvar))
  if(identical(labtext[2], NA)) ylab <- yvar else ylab <- substitute(y ~ yparen, list(yparen=yparen, y=yvar))

  # Initialize plot: add a 0 to make sure we can see the axis
  # Prevent NA values from influencing the scale of the plot 20200103
  ina <- is.na(X_d) | is.na(Y_d)
  if(!add) plot(type="n", c(X_d[!ina], 0), c(Y_d[!ina], 0), xlab=xlab, ylab=ylab, ...)
  # Add a reference rectangle
  if(plot.rect) rect(-0.01, -0.01, 0.01, 0.01, col="grey80", lwd=0)
  # show axis lines
  abline(h=0, lty=3, col = "gray30")
  abline(v=0, lty=3, col = "gray30")
  if(oldstyle) {
    # Show drop lines: dotted/solid if p-value/effect size meet criteria
    lty.X <- ifelse(abs(X_e - 50) >= 10, 1, ifelse(X_p < 0.05, 2, 0))
    lty.Y <- ifelse(abs(Y_e - 50) >= 10, 1, ifelse(Y_p < 0.05, 2, 0))
    for(i in seq_along(X_d)) {
      lines(rep(X_d[i], 2), c(0, Y_d[i]), lty=lty.X[i])
      lines(c(0, X_d[i]), rep(Y_d[i], 2), lty=lty.Y[i])
    } 
    # Point symbols: open circle, filled circle, filled square (0, 1 or 2 vars with p-value < 0.05)
    p_signif <- rowSums(data.frame(X_p < 0.05, Y_p < 0.05))
    pch <- ifelse(p_signif==2, 15, ifelse(p_signif==1, 19, 21))
  }

  # Plot points with specified color (and point style, only for oldstyle = FALSE)
  col <- rep(col, length.out=nrow(comptab))
  pch <- rep(pch, length.out=nrow(comptab))
  points(X_d, Y_d, pch=pch, col=col, bg="white", cex=cex)
  if(!identical(pt.text, NA) | !identical(pt.text, FALSE)) {
    # Add letters; for dark background, use white letters
    col[pch %in% c(15, 19)] <- "white"
    text(X_d, Y_d, pt.text[seq_along(X_d)], col=col, cex=cex.text)
  }

  # Contour 2-D kernel density estimate 20190329
  # https://stats.stackexchange.com/questions/31726/scatterplot-with-contour-heat-overlay
  if(!oldstyle & any(contour)) {
    # Include points specified by 'contour' 20191102
    Xcont <- X_d[contour]
    Ycont <- Y_d[contour]
    # Remove NA points (possible with PS (phylostrata)) 20191127
    iNA <- is.na(Xcont) | is.na(Ycont)
    Xcont <- Xcont[!iNA]
    Ycont <- Ycont[!iNA]
    if(length(Xcont) > 0) {
      dens <- kde2d(Xcont, Ycont, n = 200)
      # Add contour around 50% of points (or other fractions specified by 'probs') 20191126
      # https://stackoverflow.com/questions/16225530/contours-of-percentiles-on-level-plot
      # (snippet from emdbook::HPDregionplot from @benbolker)
      dx <- diff(dens$x[1:2])
      dy <- diff(dens$y[1:2])
      sz <- sort(dens$z)
      c1 <- cumsum(sz) * dx * dy
      levels <- sapply(probs, function(x) {
        approx(c1, sz, xout = 1 - x)$y
      })
      # Use lty = 2 and lwd = 1 if points are being plotted, or lty = 1 and lwd = 2 otherwise 20191126
      if(identical(pch[1], NA)) lty <- 1 else lty <- 2
      if(identical(pch[1], NA)) lwd <- 2 else lwd <- 1
      # Don't try to plot contours for NA levels 20191207
      if(!any(is.na(levels))) contour(dens, drawlabels = FALSE, levels = levels, lty = lty, lwd = lwd, add = TRUE, col = col.contour)
    }
  }

}

# Text for figure labels
# Moved from internal.R and exported 20200204
cplab <- list(
  nH2O = expression(italic(n)[H[2]*O]),
  DnH2O = expression(Delta*italic(n)[H[2]*O]),
  nO2 = expression(italic(n)[O[2]]),
  DnO2 = expression(Delta*italic(n)[O[2]]),
  ZC = expression(italic(Z)[C]),
  DZC = expression(Delta*italic(Z)[C]),
  logfO2 = expression(log~italic("f")[O[2]]),
  logaH2O = expression(log~italic("a")[H[2]*O]),
  nC = expression(italic(n)[C] * "/AA"),
  nN = expression(italic(n)[N] * "/AA"),
  nS = expression(italic(n)[S] * "/AA"),
  DnC = expression(Delta*italic(n)[C] * "/AA"),
  DnN = expression(Delta*italic(n)[N] * "/AA"),
  DnS = expression(Delta*italic(n)[S] * "/AA"),
  V0 = expression(list(italic("V") * degree, "cm" ^ 3 ~ "mol" ^ -1)),
  DV0 = expression(list(Delta * italic("V") * degree, "cm" ^ 3 ~ "mol" ^ -1)),
  nAA = expression(italic(n)[AA]),
  DnAA = expression(Delta*italic(n)[AA]),
  GRAVY = "GRAVY",
  DGRAVY = expression(Delta*"GRAVY"),
  pI = "pI",
  DpI = expression(Delta*"pI"),
  MW = expression("MW"),
  DMW = expression(Delta * "MW")
)
