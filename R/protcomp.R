# canprot/R/protcomp.R
# function to read protein data and calculate compositional parameters
# 20160705 jmd

protcomp <- function(uniprot = NULL, basis = getOption("basis"), aa = NULL, aa_file = NULL) {
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
  # protein formula, average oxidation state of carbon
  protein.formula <- CHNOSZ::protein.formula(aa)
  ZC <- CHNOSZ::ZC(protein.formula)
  # basis species for proteins, protein length, basis species in residue
  if(basis=="QEC") basis(c("glutamine", "glutamic acid", "cysteine", "H2O", "oxygen")) 
  else if(basis=="MTa") basis(c("methionine", "threonine", "acetic acid", "H2O", "oxygen")) 
  else if(basis=="CRa") basis(c("cysteine", "arginine", "acetic acid", "H2O", "oxygen")) 
  else basis(basis)
  protein.basis <- protein.basis(aa)
  protein.length <- protein.length(aa)
  residue.basis <- protein.basis / protein.length
  # residue formula
  residue.formula <- protein.formula / protein.length
  # return data
  out <- list(uniprot = uniprot, protein.formula=protein.formula, ZC=ZC, protein.basis=protein.basis,
    protein.length=protein.length, residue.basis=residue.basis, residue.formula=residue.formula, aa=aa)
  return(out)
}
