# canprot/R/check_IDs.R
# function to identify known UniProt IDs
# 20160703 jmd

check_IDs <- function(dat, IDcol, aa_file = NULL, updates_file = NULL) {
  # the input IDs that are NA
  input.NA <- is.na(dat[, IDcol]) | dat[, IDcol] == ""
  # the candidate IDs separated into a list
  ID_list <- strsplit(dat[, IDcol], ";")
  # the list of IDs as a vector
  ID <- unlist(ID_list)
  # get the UniProt ID in case we have e.g. sp|P62308|RUXG_HUMAN
  # for NJVS19 dataset 20191226
  if(any(grepl("\\|", ID))) ID <- sapply(strsplit(ID, "\\|"), "[", 2)
  # human proteins
  aa <- get("human_aa", human)
  # add amino acid compositions from external file if specified
  if(!is.null(aa_file)) {
    aa_dat <- read.csv(aa_file, as.is=TRUE)
    aa <- rbind(aa_dat, aa)
  }
  # assemble known IDs
  knownIDs <- sapply(strsplit(aa$protein, "|", fixed = TRUE), "[", 2)
  # if that is NA (i.e. no | separator is present) use the entire string
  ina <- is.na(knownIDs)
  knownIDs[ina] <- aa$protein[ina]
  # also include obsolete UniProt ID
  updates <- get("uniprot_updates", human)
  if(!is.null(updates_file)) {
    updates_dat <- read.csv(updates_file, as.is = TRUE)
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
  # the output IDs that are NA
  output.NA <- is.na(ID) | ID == ""
  # print which IDs became NA 20191127
  if(sum(output.NA) > sum(input.NA)) {
    new.NA <- output.NA & !input.NA
    NA.IDs <- dat[new.NA, IDcol]
    print(paste("check_IDs:", sum(new.NA), "unavailable UniProt IDs:", paste(NA.IDs, collapse = " ")))
  }
  # now apply the updates 20191119
  iold <- match(ID, updates$old)
  if(any(!is.na(iold))) {
    oldIDs <- updates$old[na.omit(iold)]
    newIDs <- updates$new[na.omit(iold)]
    print(paste("check_IDs: updating", sum(!is.na(iold)), "old UniProt IDs:", paste(oldIDs, collapse = " ")))
    ID[!is.na(iold)] <- newIDs
  }
  dat[, IDcol] <- ID
  dat
}
