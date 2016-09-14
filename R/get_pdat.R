# canprot/R/get_pdat.R
# get protein expression and composition data
# 20160712 jmd

get_pdat <- function(dataset=NULL, pdat_fun="pdat_CRC", basis="AA") {
  # datasets available from all functions
  datasets <- lapply(pdat_fun, function(x) get(x)())
  names(datasets) <- pdat_fun
  if(is.null(dataset)) return(datasets)
  # which function (if any) provides this dataset
  matchdat <- lapply(datasets, function(x) match(dataset, x))
  hasdat <- sapply(matchdat, function(x) !is.na(x))
  if(sum(hasdat) < 1) stop(paste("no function provides data for", dataset))
  if(sum(hasdat) > 1) stop(paste0("multiple functions (", paste(pdat_fun[hasdat], collapse=", "), ") provide data for ", dataset))
  # in case the name of the dataset has an =XX suffix, remove it here
  dataset <- strsplit(dataset, "=")[[1]][1]
  return(get(pdat_fun[hasdat])(dataset, basis))
}
