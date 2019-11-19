# canprot/R/check_IDs.R
# function to identify known UniProt IDs
# 20160703 jmd

check_IDs <- function(dat, IDcol, aa_file=NULL, updates_file=NULL) {
  # the candidate IDs separated into a list
  ID_list <- strsplit(dat[, IDcol], ";")
  # the list of IDs as a vector
  ID <- unlist(ID_list)
  # human proteins
  aa <- get("human_aa", canprot)
  # add amino acid compositions from external file if specified
  if(!is.null(aa_file)) {
    aa_dat <- read.csv(aa_file, as.is=TRUE)
    aa <- rbind(aa_dat, aa)
  }
  # find known IDs
  knownIDs <- sapply(strsplit(aa$protein, "|", fixed=TRUE), "[", 2)
  # if that is NA (i.e. no | separator is present) use the entire string
  ina <- is.na(knownIDs)
  knownIDs[ina] <- aa$protein[ina]
  # also check old to new UniProt ID mapping
  updates <- get("uniprot_updates", canprot)
  if(!is.null(updates_file)) {
    updates_dat <- read.csv(updates_file, as.is=TRUE)
    updates <- rbind(updates_dat, updates)
  }
  knownIDs <- c(knownIDs, updates$old)
  # check if the candidate IDs are known
  known <- match(ID, knownIDs)
  known_IDs <- ID[known > 0]
  # get the IDs back into list form
  known_IDs <- relist(known_IDs, ID_list)
  # take the first (non-NA) match 20191119
  ID <- sapply(sapply(known_IDs, na.omit), "[", 1)
  dat[, IDcol] <- ID
  dat
}
