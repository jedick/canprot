# canprot/R/remove_entries.R
# utility function to remove selected entries and print a message
# 20160713 jmd

remove_entries <- function(dat, irm, dataset, description) {
  if(is.null(dataset)) dataset <- "remove_entries"
  if(sum(irm) > 0) {
    dat <- dat[!irm, ]
    print(paste0(dataset, ": dropping ", sum(irm), " ", description, " proteins"))
  }
  return(dat)
}
