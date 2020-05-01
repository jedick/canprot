# canprot/metrics.R
# calculate various metrics from amino acid composition of proteins
# 20191027

# calculate carbon oxidation state for amino acid compositions 20180228
ZCAA <- function(AAcomp, nothing=NULL) {
  # a dummy second argument is needed because of how this function is used in JMDplots::plotMG
  # the number of carbons of the amino acids
  nC_AA <- c(Ala = 3, Cys = 3, Asp = 4, Glu = 5, Phe = 9, Gly = 2, His = 6, 
    Ile = 6, Lys = 6, Leu = 6, Met = 5, Asn = 4, Pro = 5, Gln = 5, 
    Arg = 6, Ser = 3, Thr = 4, Val = 5, Trp = 11, Tyr = 9)
  # the Ztot of the amino acids == CHNOSZ::ZC(info(info(aminoacids("")))$formula) * nC_AA
  Ztot_AA <- c(Ala = 0, Cys = 2, Asp = 4, Glu = 2, Phe = -4, Gly = 2, His = 4, 
    Ile = -6, Lys = -4, Leu = -6, Met = -2, Asn = 4, Pro = -2, Gln = 2, 
    Arg = 2, Ser = 2, Thr = 0, Val = -4, Trp = -2, Tyr = -2)
  # the ZC of the amino acids == CHNOSZ::ZC(info(info(aminoacids("")))$formula)
  ZC_AA <- Ztot_AA / nC_AA
  # find columns with names for the amino acids
  isAA <- colnames(AAcomp) %in% c("Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", 
    "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr")
  iAA <- match(colnames(AAcomp)[isAA], names(ZC_AA))
  # calculate the nC for all occurrences of each amino acid
  multC <- t(t(AAcomp[, isAA]) * nC_AA[iAA])
  # multiply nC by ZC
  multZC <- t(t(multC) * ZC_AA[iAA])
  # calculate the total ZC and nC, then the overall ZC
  ZCtot <- rowSums(multZC)
  nCtot <- rowSums(multC)
  ZCtot / nCtot
}

# calculate stoichiometric water content for amino acid compositions 20181228
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
  if(basis == "rQEC") {
    # round(residuals(lm(nH2O_AA ~ ZC(species()$ispecies))), 3)
    nH2O_AA <- c(Ala = 0.724, Cys = 0.33, Asp = 0.233, Glu = 0.248, Phe = -2.213,
      Gly = 0.833, His = -1.47, Ile = 1.015, Lys = 1.118, Leu = 1.015,
      Met = 0.401, Asn = 0.233, Pro = 0.001, Gln = 0.248, Arg = 0.427,
      Ser = 0.93, Thr = 0.924, Val = 0.877, Trp = -3.732, Tyr = -2.144)
    # subtract a constant to make the mean for human proteins = 0 20191114
    nH2O_AA <- nH2O_AA - 0.355
  }
  # water stoichiometry of amino acid biosynthetic reactions from metabolic precursors 20191205
  if(basis == "biosynth") {
    nH2O_AA <- c(Ala = 0, Cys = -1, Asp = 0, Glu = 0, Phe = -1, Gly = -1, His = -4,
      Ile = 3, Lys = 2, Leu = 3, Met = 1, Asn = -1, Pro = 0, Gln = -1,
      Arg = -2, Ser = 0, Thr = 1, Val = 2, Trp = -2, Tyr = -1)
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

# get stoichiometric O2 coefficients in amino acid biosynthetic reactions from metabolic precursors 20191205
O2AA <- function(AAcomp, basis = "biosynth") {
  if(basis == "biosynth") {
    nO2_AA <- c(Ala = -0.5, Cys = 0, Asp = -0.5, Glu = -0.5, Phe = -0.5, Gly = 1, 
      His = 0, Ile = -5, Lys = -4.5, Leu = -5, Met = -3, Asn = -0.5, 
      Pro = -1.5, Gln = -0.5, Arg = -1.5, Ser = 0, Thr = -1.5, Val = -3.5, 
      Trp = -2, Tyr = 0)
  }
  # find columns with names for the amino acids
  isAA <- colnames(AAcomp) %in% names(nO2_AA)
  iAA <- match(colnames(AAcomp)[isAA], names(nO2_AA))
  # calculate total number of O2 in reactions to form proteins
  nO2 <- rowSums(t(t(AAcomp[, isAA]) * nO2_AA[iAA]))
  # divide by number of residues (length of protein)
  nO2 / rowSums(AAcomp[, isAA])
}

# calculate GRAVY for amino acid compositions 20191024
GRAVY <- function(AAcomp) {
  # values of the hydropathy index from Kyte and Doolittle, 1982
  # doi:10.1016/0022-2836(82)90515-0
  Hind <- c(Ala =  1.8, Cys =  2.5, Asp = -3.5, Glu = -3.5, Phe =  2.8,
            Gly = -0.4, His = -3.2, Ile =  4.5, Lys = -3.9, Leu =  3.8,
            Met =  1.9, Asn = -3.5, Pro = -1.6, Gln = -3.5, Arg = -4.5,
            Ser = -0.8, Thr = -0.7, Val =  4.2, Trp = -0.9, Tyr = -1.3)
  # find columns with names for the amino acids
  isAA <- colnames(AAcomp) %in% names(Hind)
  iAA <- match(colnames(AAcomp)[isAA], names(Hind))
  # calculate total of hydropathy values for each protein
  sumHind <- rowSums(t(t(AAcomp[, isAA]) * Hind[iAA]))
  # divide by length of proteins to get grand average of hydropathy (GRAVY)
  sumHind / rowSums(AAcomp[, isAA])
}

# calculate isoelectric point for proteins 20191026
pI <- function(AAcomp) {
  # a function to calculate isoelectric point for a single amino acid composition
  onepI <- function(AA) {
    # find the column names of AAcomp that are in Ztab
    isZ <- names(AA) %in% dimnames(Ztab)[[2]]
    iZ <- match(names(AA)[isZ], dimnames(Ztab)[[2]])
    # calculate the total charge as a function of pH
    # ... the "else" is in case we have a data frame (used when first writing this function)
    if(is.numeric(AA)) Ztot <- Ztab[, iZ] %*% AA[isZ]
    else Ztot <- Ztab[, iZ] %*% as.matrix(t(AA[, isZ]))
    # find pH where charge is closest to zero
    # (absolute charge is minimized)
    ipH <- which.min(abs(Ztot))
    Ztab[ipH, 1]
  }
  # number of N- and C-terminal groups is 1, unless the input data frame has a value for number of chains
  Nterm <- Cterm <- 1
  if(!is.null(AAcomp$chains)) Nterm <- Cterm <- AAcomp$chains
  if(!"Nterm" %in% names(AAcomp)) AAcomp <- cbind(AAcomp, Nterm = Nterm)
  if(!"Cterm" %in% names(AAcomp)) AAcomp <- cbind(AAcomp, Cterm = Cterm)
  # NOTE: apply() converts the input to matrix,
  # so we extract the numeric columns of AAcomp to avoid possible coercion of all values to character
  isnum <- unlist(lapply(AAcomp, "class")) %in% c("integer", "numeric")
  myAA <- AAcomp[, isnum]
  # run the calculation for each composition
  apply(myAA, 1, onepI)
}

# calculate average molecular weight per amino acid 20200501
MWAA <- function(AAcomp) {
  # mass per residue:
  # MW_AA <- sapply(CHNOSZ::makeup(info(aminoacids(""))), mass) - mass("H2O")
  # names(MW_AA) <- aminoacids(3)
  MW_AA <- c(Ala = 71.0788, Cys = 103.1388, Asp = 115.0886, Glu = 129.11548, 
    Phe = 147.17656, Gly = 57.05192, His = 137.14108, Ile = 113.15944, 
    Lys = 128.17408, Leu = 113.15944, Met = 131.19256, Asn = 114.10384, 
    Pro = 97.11668, Gln = 128.13072, Arg = 156.18748, Ser = 87.0782, 
    Thr = 101.10508, Val = 99.13256, Trp = 186.2132, Tyr = 163.17596
  )
  # find columns with names for the amino acids
  isAA <- colnames(AAcomp) %in% names(MW_AA)
  iAA <- match(colnames(AAcomp)[isAA], names(MW_AA))
  # calculate total MW of residues in each protein
  MW <- rowSums(t(t(AAcomp[, isAA]) * MW_AA[iAA]))
  # add terminal H2O
  MW <- MW + 18.01528
  # divide by number of residues (length of protein)
  MW / rowSums(AAcomp[, isAA])
}

#########################
### UNEXPORTED OBJECT ###
### ( used in pI() )  ###
#########################

# tabulate charges for sidechains and terminal groups from pH 0 to 14
Ztab <- local({
  # a function to calculate charge as a function of pH for a single group
  ZpH <- function(pK, Z, pH) {
    alpha <- 1/(1 + 10^(Z * (pH - pK)))
    alpha * Z
  }
  # list the pKs of the groups
  pK <- list(Cterm = 3.55, Nterm = 7.5,
    Asp = 4.05, Glu = 4.45, His = 5.98,
    Cys = 9, Lys = 10, Tyr = 10, Arg = 12
  )
  # list the unit charges of the groups
  Z <- list(Cterm = -1, Nterm = 1,
    Asp = -1, Glu = -1, His = 1,
    Cys = -1, Lys = 1, Tyr = -1, Arg = 1
  )
  # get the charges for a range of pH values
  pH <- seq(0, 14, 0.01)
  Ztab <- mapply(ZpH, pK = pK, Z = Z, MoreArgs = list(pH = pH))
  # add a column with the pH values
  cbind(pH = pH, Ztab)
})
