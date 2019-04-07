# canprot/R/lapply_canprot.R
# conditionally run computations in parallel
# 20160715 jmd

lapply_canprot <- function(X, FUN, ..., varlist=NULL, min.length=10) {
  # set up worker environments and run parLapply if
  # length(X) >= min.length; otherwise, run lapply
  if(length(X) >= min.length) {
    # Use option mc.cores to choose an appropriate cluster size.
    # and set max at 2 for now (per CRAN policies)
    nCores <- min(getOption("mc.cores"), 2)
    # don't load methods package, to save startup time - ?makeCluster
    cl <- parallel::makeCluster(nCores, methods=FALSE)
    # load CHNOSZ
    message("lapply_canprot: loading CHNOSZ")
    parallel::clusterEvalQ(cl, library("CHNOSZ"))
    # load canprot
    message("lapply_canprot: loading canprot")
    parallel::clusterEvalQ(cl, library("canprot"))
    if(!is.null(varlist)) parallel::clusterExport(cl, varlist)
    # run the calculations
    message(paste("lapply_canprot: running", length(X), "calculations on", nCores, "cores"))
    out <- parallel::parLapply(cl, X, FUN, ...)
    parallel::stopCluster(cl)
  } else out <- lapply(X, FUN, ...)
  return(out)
}
