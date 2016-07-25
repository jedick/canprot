# canprot/R/protcomp.R
# function to read protein data and calculate compositional parameters
# 20160705 jmd

protcomp <- function(uniprot=NULL, ip=NULL, basis="AA", aa_file=NULL, updates_file=NULL) {
  if(is.null(ip)) {
    # get amino acid compositions of human proteins
    aa <- get("human_aa", "canprot")
    # add amino acid compositions from external file if specified
    if(!is.null(aa_file)) {
      aa_dat <- read.csv(aa_file, as.is=TRUE)
      print(paste("protcomp: adding", nrow(aa_dat), "proteins from", aa_file))
      aa <- rbind(aa_dat, aa)
    }
    if(is.null(uniprot)) {
      stop("'uniprot' is NULL")
    } else {
      # convert old to new uniprot IDs
      updates <- get("uniprot_updates", "canprot")
      # include updates from external file if specified
      if(!is.null(updates_file)) {
        updates_dat <- read.csv(updates_file, as.is=TRUE)
        updates <- rbind(updates_dat, updates)
      }
      iold <- match(uniprot, updates$old)
      if(any(!is.na(iold))) {
        oldIDs <- updates$old[na.omit(iold)]
        newIDs <- updates$new[na.omit(iold)]
        print(paste("protcomp: updating", sum(!is.na(iold)), "old UniProt IDs:", paste(oldIDs, collapse=" ")))
        # check if new IDs are duplicated
        idup <- newIDs %in% uniprot
        if(any(idup)) print(paste("protcomp: new IDs in updates",
          paste(oldIDs[idup], newIDs[idup], sep="->", collapse=" "), "are duplicated in dataset"))
        uniprot[!is.na(iold)] <- newIDs
      }
      # find the proteins listed in 'uniprot'
      alluni <- sapply(strsplit(aa$protein, "|", fixed=TRUE), "[", 2)
      iuni <- match(uniprot, alluni)
      # stop with error if any are not found
      if(any(is.na(iuni))) stop(paste("uniprot IDs not found:", paste(uniprot[is.na(iuni)], collapse=" ")))
      # warn if any are duplicated
      if(any(duplicated(iuni))) warning(paste("some uniprot IDs are duplicated:",
        paste(uniprot[duplicated(iuni)], collapse=" ")), immediate.=TRUE)
      aa <- aa[iuni, ]
    }
  } else {
    # ip is given, for pre-loaded proteins, used by canstab()
    aa <- ip
  }
  # protein formula, average oxidation state of carbon
  protein.formula <- protein.formula(aa)
  ZC <- ZC(protein.formula)
  # basis species in protein, protein length, basis species in residue
  # TODO: use CHNOSZ::basis instead (requires update to code in CHNOSZ)
  setbasis(basis)
  protein.basis <- protein.basis(aa)
  protein.length <- protein.length(aa)
  residue.basis <- protein.basis / protein.length
  # residue formula
  residue.formula <- protein.formula / protein.length
  # return data
  out <- list(protein.formula=protein.formula, ZC=ZC, protein.basis=protein.basis,
    protein.length=protein.length, residue.basis=residue.basis, residue.formula=residue.formula, aa=aa)
  return(out)
}

