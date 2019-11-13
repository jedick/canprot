# function to update UniProt IDs
# extracted from protcomp() 20191113

update_IDs <- function(dat, IDcol, updates_file = NULL) {
  # convert old to new uniprot IDs
  updates <- get("uniprot_updates", canprot)
  # include updates from external file if specified
  if(!is.null(updates_file)) {
    updates_dat <- read.csv(updates_file, as.is=TRUE)
    updates <- rbind(updates_dat, updates)
  }
  iold <- match(dat[, IDcol], updates$old)
  if(any(!is.na(iold))) {
    oldIDs <- updates$old[na.omit(iold)]
    newIDs <- updates$new[na.omit(iold)]
    print(paste("update_IDs: updating", sum(!is.na(iold)), "old UniProt IDs:", paste(oldIDs, collapse=" ")))
#    # check if new IDs are duplicated
#    idup <- newIDs %in% dat[, IDcol]
#    if(any(idup)) print(paste("protcomp: new IDs in updates",
#      paste(oldIDs[idup], newIDs[idup], sep="->", collapse=" "), "are duplicated in dataset"))
    dat[!is.na(iold), IDcol] <- newIDs
  }
  dat
}
