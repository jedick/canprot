# canprot/R/wrdiff.R
# calculate weighted difference of sums of ranks
# 20160710 jmd

rankdiff <- function(rank1, rank2, n1=NULL, n2=NULL, as.fraction=FALSE) {
  if(!is.null(n1)) sum1 <- rank1
  else {
    sum1 <- sum(rank1)
    n1 <- length(rank1)
  }
  if(!is.null(n2)) sum2 <- rank2
  else {
    sum2 <- sum(rank2)
    n2 <- length(rank2)
  }
  # scaling to account for different numbers of proteins
  # up-expressed in normal (n1) and cancer (n2) tissue
  sum1 <- 2 * sum1 * n2 / (n1 + n2)
  sum2 <- 2 * sum2 * n1 / (n1 + n2)
  sumdiff <- sum2 - sum1
  if(!as.fraction) return(sumdiff)
  else{
    # what is the maximum possible rank difference? e.g. 20 for H.H.H.H.H.C.C.C.C
    # use numeric values to avoid maxC * n1 giving NAs
    # produced by integer overflow for large data sets (e.g. WDO+15)
    min1 <- as.numeric(sum(1:n1))
    maxC <- as.numeric(sum((1:n2) + n1))
    maxdiff <- 2 * (maxC * n1 - min1 * n2) / (n1 + n2)
    # express the rank difference as a fraction of maximum
    rankdiff <- sumdiff / maxdiff
    return(rankdiff)
  }
}
