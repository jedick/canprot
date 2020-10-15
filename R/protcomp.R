# canprot/R/protcomp.R
# function to get amino acid compositions for list of UniProt IDs
# 20160705 jmd

protcomp <- function(uniprot = NULL, aa = NULL, aa_file = NULL) {
  if(is.null(aa)) {
    # get amino acid compositions of human proteins
    aa <- get("human_aa", human)
    # add amino acid compositions from external file if specified
    if(!is.null(aa_file)) {
      aa_dat <- read.csv(aa_file, as.is=TRUE)
      print(paste("protcomp: adding", nrow(aa_dat), "proteins from", aa_file))
      aa <- rbind(aa_dat, aa)
    }
    if(is.null(uniprot)) {
      stop("'uniprot' is NULL")
    } else {
      # find the proteins listed in 'uniprot' - first look at the ID after the | separator
      alluni <- sapply(strsplit(aa$protein, "|", fixed=TRUE), "[", 2)
      # if that is NA (i.e. no | separator is present) use the entire string
      ina <- is.na(alluni)
      alluni[ina] <- aa$protein[ina]
      iuni <- match(uniprot, alluni)
      # stop with error if any are not found
      if(any(is.na(iuni))) stop(paste("uniprot IDs not found:", paste(uniprot[is.na(iuni)], collapse=" ")))
      # warn if any are duplicated
      if(any(duplicated(iuni))) warning(paste("some uniprot IDs are duplicated:",
        paste(uniprot[duplicated(iuni)], collapse=" ")), immediate.=TRUE)
      aa <- aa[iuni, ]
    }
  }
  # return data
  out <- list(uniprot = uniprot, aa = aa)
  return(out)
}
