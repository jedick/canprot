# canprot/qdiff.R
# 20200505 jmd

# quantile plot for up- and down-regulated proteins
# first version 20200428
qdiff <- function(pd = pdat_colorectal("JKMF10"), vars = c("ZC", "nH2O")) {
  # initialize plot
  if(length(vars)==2) par(mfrow = c(2, 1))
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
    plot(range(up, dn), c(0, 1), type = "n", xlab = xlab, ylab = "quantile")
    # plot quantiles
    lines(sort(up), (1:length(up)) / (length(up)+1), col = "red", lty = 2)
    lines(sort(dn), (1:length(dn)) / (length(dn)+1))
    # draw a line at the 0.5 quantile
    lines(c(median(dn), median(up)), c(0.5, 0.5))
  }
}
