# canprot/R/cleanup.R
# Function to clean up data (remove proteins with unavailable IDs,
# ambiguous expression ratios, and duplicated IDs)
# 20190407 jmd

cleanup <- function(dat, IDcol, up2 = NULL) {

  # Utility function to remove selected entries and print a message 20160713
  remove_entries <- function(dat, irm, description) {
    if(sum(irm) > 0) {
      dat <- dat[!irm, , drop = FALSE]
      if(sum(irm)==1) ptxt <- "protein" else ptxt <- "proteins"
      print(paste("cleanup: removing", sum(irm), description, ptxt))
    }
    return(dat)
  }

  if(is.logical(IDcol)) {
    dat <- remove_entries(dat, IDcol, "selected")
  } else {
    # Add up2 column to dat
    if(!is.null(up2)) dat <- cbind(dat, up2 = up2)
    # Which column has the IDs
    if(is.character(IDcol)) IDcol <- match(IDcol, colnames(dat))
    # Drop NA or "" IDs
    unav <- is.na(dat[, IDcol]) | dat[, IDcol] == ""
    dat <- remove_entries(dat, unav, "unavailable")
    if(!is.null(up2)) {
      # Drop unquantified proteins 20191120
      dat <- remove_entries(dat, is.na(dat$up2), "unquantified")
      # Drop proteins with ambiguous expression ratios
      up <- dat[dat$up2, IDcol]
      down <- dat[!dat$up2, IDcol]
      iambi <- dat[, IDcol] %in% intersect(up, down)
      dat <- remove_entries(dat, iambi, "ambiguous")
    }
    # Drop duplicated proteins
    dat <- remove_entries(dat, duplicated(dat[, IDcol]), "duplicated")
  }
  dat
}
