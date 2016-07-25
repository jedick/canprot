# canprot/R/rankplot.R
# make logaH2O - logfO2 diagrams of ranking of chemical affinity
# 20160710 jmd

rankplot <- function(pdat, T=37, plot.rank=TRUE, main=NULL, res=300) {
  # get protein data
  aa <- pdat$pcomp$aa
  ip <- add.protein(aa)
  # assemble limits of logfO2, logaH2O
  # use res+1 here to ensure that x- and y-directions aren't accidentally transposed
  # (would cause an error in image())
  H2O <- c(-10, 10, res)
  if(pdat$basis=="inorganic") O2 <- c(-80, -60, res+1)
  else O2 <- c(-75, -55, res+1)
  ys <- seq(H2O[1], H2O[2], length.out=H2O[3])
  xs <- seq(O2[1], O2[2], length.out=O2[3])
  # what is the log activity of protein corresponding to unit activity of residues?
  loga.protein <- log10(1/pdat$pcomp$protein.length)
  # calculate affinities
  a <- affinity(O2=O2, H2O=H2O, T=T, iprotein=ip, loga.protein=loga.protein)
  if(!plot.rank) {
    col <- ifelse(pdat$up2, cpcol$red, cpcol$blue)
    d <- diagram(a, names=pdat$names, fill=col, as.residue=TRUE, tplot=FALSE)
    # redraw box because it gets obscured by the fill
    box()
  } else {
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
    rankdiff <- 100*rankdiff(rank_H, rank_C, n_H, n_C, as.fraction=TRUE)
    print(paste("weighted rank difference % range", paste(range(rankdiff), collapse=" ")))
    # display greater healthy and cancer ranks by colored zones on image
    dcol <- colorspace::diverge_hcl(1000, c = 100, l = c(50, 90), power = 1)
    if(any(rankdiff < 0)) {
      # so that anything over 50% is deep blue
      blues <- rev(c(rep(dcol[1], 500), dcol[1:500]))
      # select range of colors corresponding to rank difference percentage
      blues <- blues[1:-round(10*range(rankdiff)[1])]
    } else blues <- character()
    if(any(rankdiff > 0)) {
      # so that anything over 50% is deep red
      reds <- c(dcol[501:1000], rep(dcol[1000], 500))
      reds <- reds[1:round(10*range(rankdiff)[2])]
    } else reds <- character()
    col <- c(rev(blues), reds)
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
