# canprot/R/diffplot.R
# plot differences in means of ZC and nH2O
# and indicate low p-values with lines
# 20160715 jmd

diffplot <- function(comptab, col="black", plot.rect=FALSE) {
  # convert to data frame if needed
  if(!is.data.frame(comptab)) comptab <- do.call(rbind, comptab)
  # get mean difference, common language effect size and p-value
  ZC_d <- comptab$ZC.diff
  ZC_e <- signif(comptab$ZC.CLES, 2)
  ZC_p <- comptab$ZC.p.value
  nH2O_d <- comptab$nH2O.diff
  nH2O_e <- signif(comptab$nH2O.CLES, 2)
  nH2O_p <- comptab$nH2O.p.value
  # set up plot
  xlab <- substitute("mean difference (" * x * ")", list(x=cplab$ZC[[1]]))
  ylab <- substitute("mean difference (" * y * ")", list(y=cplab$nH2O[[1]]))
  # initialize plot: add a 0 to make sure we can see the axis
  plot(type="n", c(ZC_d, 0), c(nH2O_d, 0), xlab=xlab, ylab=ylab)
  # add a reference rectangle
  if(plot.rect) rect(-0.01, -0.01, 0.01, 0.01, col="grey80", lwd=0)
  # show axis lines
  abline(h=0, lty=3)
  abline(v=0, lty=3)
  # show drop lines: dotted/solid if p-value/effect size meet criteria
  lty.ZC <- ifelse(abs(ZC_e - 50) >= 10, 1, ifelse(ZC_p < 0.05, 2, 0))
  lty.nH2O <- ifelse(abs(nH2O_e - 50) >= 10, 1, ifelse(nH2O_p < 0.05, 2, 0))
  for(i in seq_along(ZC_d)) {
    lines(rep(ZC_d[i], 2), c(0, nH2O_d[i]), lty=lty.ZC[i])
    lines(c(0, ZC_d[i]), rep(nH2O_d[i], 2), lty=lty.nH2O[i])
  } 
  # point symbols: open circle, filled circle, filled square (0, 1 or 2 vars with p-value < 0.05)
  p_signif <- rowSums(data.frame(ZC_p < 0.05, nH2O_p < 0.05))
  pch <- ifelse(p_signif==2, 15, ifelse(p_signif==1, 19, 21))
  # plot points with specified color
  col <- rep(col, length.out=nrow(comptab))
  points(ZC_d, nH2O_d, pch=pch, col=col, bg="white")
}
