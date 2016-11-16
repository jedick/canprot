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
    msgout("lapply_canprot: loading CHNOSZ\n")
    parallel::clusterEvalQ(cl, library("CHNOSZ"))
    # we don't do this because it's implied by data(canprot):
    #parallel::clusterCall(cl, "data", list="thermo")
    # load canprot
    msgout("lapply_canprot: loading canprot and setting up canprot environment\n")
    parallel::clusterEvalQ(cl, library("canprot"))
    # we don't do this because a simple data(canprot) is nearly as fast
    #parallel::clusterExport(cl, c("human_aa", "uniprot_updates"), as.environment("canprot"))
    #parallel::clusterCall(cl, "attach", what=NULL, name="canprot")
    #parallel::clusterEvalQ(cl, {
    #  assign("human_aa", get("human_aa"), "canprot")
    #  assign("uniprot_updates", get("uniprot_updates"), "canprot")
    #})
    parallel::clusterCall(cl, "data", list="canprot")
    ## export the user's variables
    #if(! "" %in% varlist) {
    #  parallel::clusterExport(cl, varlist, envir)
    #  msgout(paste("lapply_canprot: exporting variable(s) from",
    #    attr(envir, "name"), "environment:", paste(varlist, collapse=", "), "\n"))
    #}
    if(!is.null(varlist)) parallel::clusterExport(cl, varlist)
    # run the calculations
    msgout(paste("lapply_canprot: running", length(X), "calculations on", nCores, "cores\n"))
    out <- parallel::parLapply(cl, X, FUN, ...)
    parallel::stopCluster(cl)
  } else out <- lapply(X, FUN, ...)
  return(out)
}
