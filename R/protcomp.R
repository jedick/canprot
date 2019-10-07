# canprot/R/protcomp.R
# function to read protein data and calculate compositional parameters
# 20160705 jmd

protcomp <- function(uniprot=NULL, ip=NULL, basis="QEC", aa_file=NULL, updates_file=NULL) {
  if(is.null(ip)) {
    # get amino acid compositions of human proteins
    aa <- get("human_aa", canprot)
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
      updates <- get("uniprot_updates", canprot)
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
  } else {
    # ip is given, for pre-loaded proteins, used by canstab()
    aa <- ip
  }
  # protein formula, average oxidation state of carbon
  protein.formula <- CHNOSZ::protein.formula(aa)
  ZC <- CHNOSZ::ZC(protein.formula)
  # basis species for proteins, protein length, basis species in residue
  if(basis=="rQEC") basis("QEC") else basis(basis)
  protein.basis <- protein.basis(aa)
  protein.length <- protein.length(aa)
  residue.basis <- protein.basis / protein.length
  if(basis=="rQEC") residue.basis[, "H2O"] <- H2OAA(aa, basis)
  # residue formula
  residue.formula <- protein.formula / protein.length
  # return data
  out <- list(protein.formula=protein.formula, ZC=ZC, protein.basis=protein.basis,
    protein.length=protein.length, residue.basis=residue.basis, residue.formula=residue.formula, aa=aa)
  return(out)
}

# calculate nH2O for amino acid compositions 20181228
# function copied from JMDplots package 20191007
H2OAA <- function(AAcomp, basis = "rQEC") {
  # how to use CHNOSZ to get the number of H2O in reactions
  # to form amino acid residues from the "QEC" basis:
  ## basis("QEC")
  ## species(aminoacids(3))
  ## nH2O_AA <- species()[["H2O"]]
  # subtract one H2O to make residues
  ## nH2O_AA <- nH2O_AA - 1
  ## names(nH2O_AA) <- aminoacids(3)
  ## dput(nH2O_AA)
  if(basis == "QEC") {
    nH2O_AA <- c( Ala = -0.4, Cys =   -1, Asp = -1.2, Glu =   -1, Phe = -3.2, Gly = -0.6, His = -2.8,
      Ile =  0.2, Lys =  0.2, Leu =  0.2, Met = -0.6, Asn = -1.2, Pro =   -1, Gln =   -1,
      Arg = -0.8, Ser = -0.4, Thr = -0.2, Val =    0, Trp = -4.8, Tyr = -3.2)
  }
  # residual water content with QEC basis
  ## round(residuals(lm(nH2O_AA ~ ZC(species()$ispecies))), 3)
  if(basis == "rQEC") {
    nH2O_AA <- c(Ala = 0.724, Cys = 0.33, Asp = 0.233, Glu = 0.248, Phe = -2.213,
      Gly = 0.833, His = -1.47, Ile = 1.015, Lys = 1.118, Leu = 1.015,
      Met = 0.401, Asn = 0.233, Pro = 0.001, Gln = 0.248, Arg = 0.427,
      Ser = 0.93, Thr = 0.924, Val = 0.877, Trp = -3.732, Tyr = -2.144)
  }
  # find columns with names for the amino acids
  isAA <- colnames(AAcomp) %in% names(nH2O_AA)
  iAA <- match(colnames(AAcomp)[isAA], names(nH2O_AA))
  # calculate total number of H2O in reactions to form proteins
  nH2O <- rowSums(t(t(AAcomp[, isAA]) * nH2O_AA[iAA]))
  # add one to account for terminal groups
  nH2O <- nH2O + 1
  # divide by number of residues (length of protein)
  nH2O / rowSums(AAcomp[, isAA])
  # to check this function:
  #  basis("QEC")
  #  H2O.ref <- protein.basis(1:6)[, "H2O"] / protein.length(1:6)
  #  AAcomp <- thermo()$protein[1:6, ]
  #  H2O.fun <- H2OAA(AAcomp, "QEC")
  #  stopifnot(H2O.ref == H2O.fun)
}

