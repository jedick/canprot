# canprot/R/rankplot.R
# make logaH2O - logfO2 diagrams of ranking of chemical affinity
# 20160710 jmd

rankplot <- function(pdat, T=37, what="rankdiff", main=NULL, res=300, plot.it=TRUE, xlim=c(-75, -55), ylim=c(-10, 10)) {
  # get protein data
  aa <- pdat$pcomp$aa
  ip <- add.protein(aa)
  # assemble limits of logfO2, logaH2O
  # use res+1 here to ensure that x- and y-directions aren't accidentally transposed
  # (would cause an error in image())
  H2O <- c(ylim, res)
  O2 <- c(xlim, res+1)
  ys <- seq(H2O[1], H2O[2], length.out=H2O[3])
  xs <- seq(O2[1], O2[2], length.out=O2[3])
  # what is the log activity of protein corresponding to unit activity of residues?
  loga.protein <- log10(1/pdat$pcomp$protein.length)
  # calculate affinities
  a <- affinity(O2=O2, H2O=H2O, T=T, iprotein=ip, loga.protein=loga.protein)
  if(identical(what, "predominance")) {
    col <- ifelse(pdat$up2, cpcol$red, cpcol$blue)
    d <- diagram(a, names=pdat$names, fill=col, as.residue=TRUE, tplot=FALSE)
    # redraw box because it gets obscured by the fill
    box()
  } else if(identical(what, "rankdiff")) {
    # plot affinity difference (healthy - cancer)
    Aarr <- list2array(a$values)
    # grid points on the plot
    grid <- expand.grid(seq_along(xs), seq_along(ys))
    rank_C <- palply("", 1:nrow(grid), function(i) {
      # the ranks at this grid point
      r <- rank(Aarr[grid[i, 1], grid[i, 2], ]/pdat$pcomp$protein.length)
      # the sum of ranks for cancer proteins at this grid point
      sum_C <- sum(r[pdat$up2])
      return(sum_C)
    })
    # numbers of cancer, healthy, and total proteins
    n_C <- sum(pdat$up2)
    n_H <- sum(!pdat$up2)
    # cancer and healthy rank sums in matrix format
    rank_C <- matrix(unlist(rank_C), nrow=length(xs), ncol=length(ys))
    rank_H <- sum(1:(n_C + n_H)) - rank_C
    # calculate weighted rank difference in percent
    rankdiff <- 100*rankdiff(rank_H, rank_C, n_H, n_C)
    print(paste("weighted rank difference % range", paste(range(rankdiff), collapse=" ")))
    if(!plot.it) return(list(xs=xs, ys=ys, rankdiff=rankdiff, xlab=cplab$logfO2, ylab=cplab$logaH2O))
    # display greater healthy and cancer ranks by colored zones on image
    col <- get_colors(rankdiff, max50=TRUE)
    image(xs, ys, rankdiff, col=col, useRaster=TRUE, xlab=cplab$logfO2, ylab=cplab$logaH2O)
    # show equal-rank line
    contour(xs, ys, rankdiff, levels=0, add=TRUE, drawlabels=FALSE)
    contour(xs, ys, rankdiff, levels=seq(-100, 100, by=10), lty=3, add=TRUE, drawlabels=FALSE)
    box()
  }
  if(is.null(main)) main <- pdat$description
  title(main=main, cex.main=1.1)
  mtext(pdat$dataset, side=4, cex=0.85, las=0, adj=0, line=-0.1)
}
