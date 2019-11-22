# canprot/R/cleanup.R
# function to clean up data (remove proteins with unavailable IDs,
# ambiguous expression ratios, and duplicated IDs)
# 20190407 jmd

cleanup <- function(dat, IDcol, up2 = NULL) {

  # 20160713 utility function to remove selected entries and print a message
  remove_entries <- function(dat, irm, description) {
    if(sum(irm) > 0) {
      dat <- dat[!irm, ]
      print(paste("cleanup: removing", sum(irm), description, "proteins"))
    }
    return(dat)
  }

  if(is.logical(IDcol)) {
    dat <- remove_entries(dat, IDcol, "selected")
  } else {
    # add up2 column to dat
    if(!is.null(up2)) dat <- cbind(dat, up2 = up2)
    # which column has the IDs
    if(is.character(IDcol)) IDcol <- match(IDcol, colnames(dat))
    # drop NA or "" IDs
    unav <- is.na(dat[, IDcol]) | dat[, IDcol] == ""
    dat <- remove_entries(dat, unav, "unavailable")
    if(!is.null(up2)) {
      # drop unquantified proteins 20191120
      dat <- remove_entries(dat, is.na(dat$up2), "unquantified")
      # drop proteins with ambiguous expression ratios
      up <- dat[dat$up2, IDcol]
      down <- dat[!dat$up2, IDcol]
      iambi <- dat[, IDcol] %in% intersect(up, down)
      dat <- remove_entries(dat, iambi, "ambiguous")
    }
    # drop duplicated proteins
    dat <- remove_entries(dat, duplicated(dat[, IDcol]), "duplicated")
  }
  dat
}
