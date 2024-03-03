# canprot/human.aa.R
# Get amino acid compositions for human proteins from UniProt IDs
# 20160705 jmd

human.aa <- function(uniprot = NULL, aa_file = NULL, stop.if.missing = FALSE, warn.if.duplicated = FALSE) {
  # Get amino acid compositions of human proteins
  aa <- get("human_aa", canprot)
  # Add amino acid compositions from external file if specified
  if(!is.null(aa_file)) {
    aa_dat <- read.csv(aa_file, as.is=TRUE)
    print(paste("human.aa: adding", nrow(aa_dat), "proteins from", aa_file))
    aa <- rbind(aa_dat, aa)
  }
  if(is.null(uniprot)) {
    stop("'uniprot' is NULL")
  } else {
    # Find the proteins listed in 'uniprot' - first look at the ID after the | separator
    alluni <- sapply(strsplit(aa$protein, "|", fixed = TRUE), "[", 2)
    # If that is NA (i.e. no | separator is present) use the entire string
    ina <- is.na(alluni)
    alluni[ina] <- aa$protein[ina]
    iuni <- match(uniprot, alluni)
    if(stop.if.missing) {
      # Stop with error if any IDs are not found
      if(any(is.na(iuni))) stop(paste("uniprot IDs not found:", paste(uniprot[is.na(iuni)], collapse = " ")))
    }
    if(warn.if.duplicated) {
      # Warn if any IDs are duplicated
      if(any(duplicated(iuni))) warning(paste("some uniprot IDs are duplicated:",
        paste(uniprot[duplicated(iuni)], collapse=" ")), immediate. = TRUE)
    }
    aa <- aa[iuni, ]
  }
  # Return amino acid compositions
  aa
}
