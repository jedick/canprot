# canprot/R/diffplot.R
# plot differences in means of ZC and nH2O
# and indicate low p-values with lines
# 20160715 jmd

diffplot <- function(comptab, col="black") {
  # convert to data frame if needed
  if(!is.data.frame(comptab)) comptab <- do.call(rbind, comptab)
  # get mean difference, common language effect size and p-value
  ZC_d <- comptab$ZC.diff
  ZC_e <- signif(comptab$ZC.CLES, 2)
  ZC_p <- comptab$ZC.p.value
  nH2O_d <- comptab$nH2O.diff
  nH2O_e <- signif(comptab$nH2O.CLES, 2)
  nH2O_p <- comptab$nH2O.p.value
  # consider thresholds for either effect size or p-value
  ZC_sig <- abs(ZC_e - 50) >= 10 | ZC_p < 0.05
  nH2O_sig <- abs(nH2O_e - 50) >= 10 | nH2O_p < 0.05
  # point symbols: open circle (0 vars), filled circle (1 var), filled square (2 vars)
  pch_open <- ifelse(!(ZC_sig | nH2O_sig), 1, NA)
  pch_circle <- ifelse(xor(ZC_sig, nH2O_sig), 19, NA)
  pch_square <- ifelse(ZC_sig & nH2O_sig, 15, NA)
  # set up plot
  xlab <- substitute("difference in means (" * x * ")", list(x=cplab$ZC[[1]]))
  ylab <- substitute("difference in means (" * y * ")", list(y=cplab$nH2O[[1]]))
  # add a 0 to make sure we can see the axis
  plot(type="n", c(ZC_d, 0), c(nH2O_d, 0), xlab=xlab, ylab=ylab)
  # show axis lines
  abline(h=0, lty=3)
  abline(v=0, lty=3)
  # show drop lines: dotted/solid if one/both of effect size and p-value meet criteria
  ZC_sigsig <- abs(ZC_e - 50) >= 10 & ZC_p < 0.05
  nH2O_sigsig <- abs(nH2O_e - 50) >= 10 & nH2O_p < 0.05
  for(i in seq_along(ZC_d)) {
    if(ZC_sigsig[i]) lty <- 1 else if(ZC_sig[i]) lty <- 2 else lty <- 0
    lines(rep(ZC_d[i], 2), c(0, nH2O_d[i]), lty=lty)
    if(nH2O_sigsig[i]) lty <- 1 else if(nH2O_sig[i]) lty <- 2 else lty <- 0
    lines(c(0, ZC_d[i]), rep(nH2O_d[i], 2), lty=lty)
  } 
  # plot points with specified color
  col <- rep(col, length.out=nrow(comptab))
  points(ZC_d, nH2O_d, pch=pch_open, col=col)
  points(ZC_d, nH2O_d, pch=pch_circle, col=col)
  points(ZC_d, nH2O_d, pch=pch_square, cex=1.1, col=col)
}
