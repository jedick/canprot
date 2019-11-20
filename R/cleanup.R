# canprot/R/cleanup.R
# function to clean up data (remove proteins with unavailable IDs,
# ambiguous expression ratios, and duplicated IDs)
# 20190407 jmd

cleanup <- function(dat, IDcol, dataset, up2) {
  # add up2 column to dat
  dat <- cbind(dat, up2 = up2)
  # which column has the IDs
  IDcol <- match(IDcol, colnames(dat))
  # drop NA or "" IDs
  unav <- is.na(dat[, IDcol]) | dat[, IDcol] == ""
  dat <- remove_entries(dat, unav, dataset, "unavailable")
  # drop unquantified proteins 20191120
  dat <- remove_entries(dat, is.na(dat$up2), dataset, "unquantified")
  # drop proteins with ambiguous expression ratios
  up <- dat[dat$up2, IDcol]
  down <- dat[!dat$up2, IDcol]
  iambi <- dat[, IDcol] %in% intersect(up, down)
  dat <- remove_entries(dat, iambi, dataset, "ambiguous")
  # drop duplicated proteins
  dat <- remove_entries(dat, duplicated(dat[, IDcol]), dataset, "duplicated")
}
