# canprot/potential.R
# plotting potential diagrams
# 20170611

# calculate rank potential diagrams for each dataset in a group
groupplots <- function(group="hypoxia_ZC_down", each100=FALSE, res=50, plot.it=TRUE) {
  idat <- get_idat(group)
  # make figures and return coordinates of equipotential lines
  if(plot.it) {
    # number of plots - 1 extra for title
    np <- sum(idat$idat) + 1
    if(np <= 28) mfrow <- c(4, 7)
    if(np <= 24) mfrow <- c(4, 6)
    if(np <= 20) mfrow <- c(4, 5)
    if(np <= 18) mfrow <- c(3, 6)
    if(np <= 16) mfrow <- c(4, 4)
    if(np <= 15) mfrow <- c(3, 5)
    if(np <= 12) mfrow <- c(3, 4)
    if(np <= 9) mfrow <- c(3, 3)
    if(np <= 8) mfrow <- c(2, 4)
    if(np <= 6) mfrow <- c(2, 3)
    if(np <= 4) mfrow <- c(2, 2)
    par(mfrow=mfrow)
    par(mar=c(3.5, 4, 2, 1.5))
    par(mgp=c(2.5, 1, 0))
    # first plot the title
    plot.new()
    opar <- par(xpd=TRUE)
    text(0.5, 0.8, idat$group, cex=2)
    text(0.5, 0.6, get_main(group), cex=2)
    par(opar)
  }
  # m=1 here to calculate lines for each dataset individually
  # expand x- and y-axes (48 log units)
  cpeach <- calcpot(what=idat$what, datasets=idat$datasets, res=res, plot.it=plot.it, each100=each100, xlim=c(-88, -40), ylim=c(-24, 24), m=1)
  return(list(group=group, cpeach=cpeach))
}

# plot merged ranked potential diagrams for groups of datasets
mergedplot <- function(gpresult, each100=FALSE, res=50) {
  idat <- get_idat(gpresult$group)
  cpresult <- calcpot(what=idat$what, datasets=idat$datasets, res=res, each100=each100)
  # the merged potential diagram
  plotit(cpresult[[1]][[1]], cpresult[[2]][[1]])
  # median and interquartile range of the individual equipotential lines
  plot_median(gpresult$cpeach, add=TRUE, transpose=grepl("ZC", gpresult$group))
  # plot the merged equipotential line again (for emphasis)
  plotit(cpresult[[1]][[1]], cpresult[[2]][[1]], add=TRUE)
}

### Internal Functions ###

get_idat <- function(group) {
  what <- strsplit(group, "_")[[1]][1]
  metric <- strsplit(group, "_")[[1]][2]
  direction <- strsplit(group, "_")[[1]][3]
  extdatadir <- system.file("extdata", package="canprot")
  file <- paste0(extdatadir, "/summary/summary_", what, ".csv")
  dat <- read.csv(file, as.is=TRUE)
  cnames <- colnames(dat)
  # the difference for this metric
  diff <- dat[, grepl(metric, cnames) & grepl("diff", cnames)]
  # the p-value and CLES for the other metric
  p.value <- dat[, !grepl(metric, cnames) & grepl("p.value", cnames)]
  CLES <- dat[, !grepl(metric, cnames) & grepl("CLES", cnames)]
  # the rows that have the difference in the specified direction
  # and the magnitude of the difference is at least 0.01
  # and for which the other metric doesn't have a large or significant change
  if(direction=="up") idat <- diff > 0.01 & !(p.value < 0.05 | abs(signif(CLES, 2) - 50) >= 10)
  if(direction=="down") idat <- diff < -0.01 & !(p.value < 0.05 | abs(signif(CLES, 2) - 50) >= 10)
  return(list(what=what, idat=idat, datasets=dat$dataset[idat]))
}

# function to plot values with colors and labels
plotit <- function(values, rankdat, add=FALSE) {
  col <- get_colors(values)
  with(rankdat, {
    if(!add) {
      image(xs, ys, values, col=col, useRaster=TRUE, xlab=xlab, ylab=ylab, axes=FALSE)
      if(diff(range(xs)) > 9) {
        axis(1, at=seq(-88, -40, by=12))
        axis(2, at=seq(-24, 24, by=12), las=1)
      } else {
        axis(1, at=seq(-70, -62, by=2))
        axis(2, at=seq(-6, 2, by=2), las=1)
      }
      box()
    }
    # show equal-rank line
    contour(xs, ys, values, levels=0, add=TRUE, drawlabels=FALSE, lwd=1.5, col="white")
  })
}

get_main <- function(group, n=NA) {
  what <- strsplit(group, "_")[[1]][1]
  metric <- strsplit(group, "_")[[1]][2]
  direction <- strsplit(group, "_")[[1]][3]
  if(direction=="up") than <- "> 0.01" else than <- "< -0.01"
  if(metric=="ZC") main <- substitute(Delta*italic(Z)[C]~than, list(than=than))
  if(metric=="H2O") main <- substitute(Delta*italic(bar(n))[H[2]*O]~than, list(than=than))
  if(!is.na(n)) {
    # make long titles to be placed on two lines
    if(what %in% c("hypoxia", "osmotic")) {
      if(what=="hypoxia") what <- "hypoxia or 3D culture"
      if(what=="osmotic") what <- "hyperosmotic stress"
      main <- c(what, substitute(main~~"("*x*")", list(what=what, main=main, x=n)))
    } else main <- substitute(what~~main~~"("*x*")", list(what=what, main=main, x=n))
  }
  return(main)
}

# plot median and 1st and 3rd quartiles of equipotential lines 20161021
plot_median <- function(eachplotres, each100=FALSE, transpose=FALSE, add=FALSE, type=7) {
  # show scatterplot of points on all lines
  if(!add) plot(eachplotres$x, eachplotres$y, pch=".")
  if(transpose) {
    ys <- seq(-24, 24, by=0.2)
    xs <- sapply(unique(eachplotres$comb), function(comb) {
      idat <- eachplotres$comb == comb
      # use loess to get smooth lines
      indlo <- loess(x~y, data=eachplotres[idat, ], span=0.2)
      spx <- predict(indlo, newdata=ys)
      # show loess lines
      if(!add) lines(spx, ys, col="gray")
      return(spx)
    })
    # calculate and plot median
    xmedian <- apply(xs, 1, median, na.rm=TRUE)
    lines(xmedian, ys, lwd=2)
    # calculate and plot 1st and 3rd quartiles
    xq1 <- apply(xs, 1, quantile, probs=0.25, na.rm=TRUE, type=type)
    xq3 <- apply(xs, 1, quantile, probs=0.75, na.rm=TRUE, type=type)
    lines(xq1, ys, lty=2)
    lines(xq3, ys, lty=2)
  } else {
    xs <- seq(-88, -40, by=0.2)
    ys <- sapply(unique(eachplotres$comb), function(comb) {
      idat <- eachplotres$comb == comb
      # use loess to get smooth lines
      indlo <- loess(y~x, data=eachplotres[idat, ], span=0.2)
      spy <- predict(indlo, newdata=xs)
      # show loess lines
      if(!add) lines(xs, spy, col="gray")
      return(spy)
    })
    # calculate and plot median
    ymedian <- apply(ys, 1, median, na.rm=TRUE)
    lines(xs, ymedian, lwd=1.5)
    # calculate and plot 1st and 3rd quartiles
    yq1 <- apply(ys, 1, quantile, probs=0.25, na.rm=TRUE, type=type)
    yq3 <- apply(ys, 1, quantile, probs=0.75, na.rm=TRUE, type=type)
    lines(xs, yq1, lty=2)
    lines(xs, yq3, lty=2)
  }
}

# merging potential diagrams for pancreatic cancer 20161012
calcpot <- function(what="pancreatic", datasets=c("LHE+04", "PCS+11_PDAC"),
                            basis="QEC+", res=50, plot.it=FALSE,
                            each100=FALSE, xlim=c(-70.2, -61.8), ylim=c(-6.2, 2.2), m=0) {
  # get the data
  pdat_fun <- paste0("pdat_", what)
  rankdat <- lapply_canprot(datasets, function(dataset) {
    pdat <- get_pdat(dataset, pdat_fun, basis)
    rankplot(pdat, res=res, plot.it=FALSE, xlim=xlim, ylim=ylim)
  }, varlist=pdat_fun, min.length=5)
  rankdiff <- lapply(rankdat, "[[", "rankdiff")
  extdatadir <- system.file("extdata", package="canprot")
  file <- paste0(extdatadir, "/summary/summary_", what, ".csv")
  dat <- read.csv(file, as.is=TRUE)
  # function to scale all values by same factor such that maximum is 100
  scale100 <- function(valuelist, each100=FALSE) {
    if(each100) lapply(valuelist, function(x) x * 100 / max(abs(range(x))))
    else {
      totalrange <- range(valuelist)
      lapply(valuelist, function(x) x * 100 / max(abs(totalrange)))
    }
  }
  # scale and display the rank differences
  rankdiff <- scale100(rankdiff, each100)
  # calculate the percentages of total values
  rabs <- lapply(rankdiff, abs)
  rsum <- sapply(rabs, sum)
  rperc <- round(100 * rsum / sum(rsum), 1)
  if(plot.it) for(i in seq_along(rankdat)) {
    plotit(rankdiff[[i]], rankdat[[i]])
    mtext(datasets[i], side=4, cex=0.85, las=0, adj=0, line=0.1)
    # show the dataset description and letter
    idat <- match(datasets[i], dat$dataset)
    label.figure(dat$set[idat], xfrac=0.25, cex=1.5)
    label.figure(paste0("     (", rperc[i], "%)"), adj=0, xfrac=0.25, cex=1.5)
    # put a circle around the letter
    opar <- par(xpd=NA)
    points(grconvertX(0.25, "nfc"), grconvertY(0.95, "nfc"), cex=3.5)
    par(opar)
  }
  if(m) {
    # merge the rank differences for all combinations of 'm' number of datasets
    combs <- combn(length(rankdat), m)
    comb <- x <- y <- numeric()
    print(paste("calcpot: calculating equipotential lines for", ncol(combs), "combinations"))
    for(i in 1:ncol(combs)) {
      merged <- Reduce("+", rankdiff[combs[, i]])
      cl <- contourLines(rankdat[[1]]$xs, rankdat[[1]]$ys, merged, levels=0)
      if(length(cl) > 0) {
        # in case the line has multiple segments
        clx <- do.call(c, lapply(cl, "[[", 2))
        cly <- do.call(c, lapply(cl, "[[", 3))
        x <- c(x, clx)
        y <- c(y, cly)
        comb <- c(comb, rep(i, length(clx)))
      }
    }
    out <- data.frame(comb=comb, x=x, y=y)
    return(out)
  } else {
    # merge the rank differences for all datasets
    merged <- Reduce("+", rankdiff) / length(rankdat)
    merged <- scale100(list(merged))
    return(invisible(list(merged=merged, rankdat=rankdat)))
  }
}
