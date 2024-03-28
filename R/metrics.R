# canprot/metrics.R
# Calculate various metrics from amino acid composition of proteins
# 20191027

# Keep a list of the objects in the current environment before we add the functions 20240304
cplist0 <- ls()

# Carbon oxidation state 20180228
# An unused second argument (...) is provided for flexible do.call() constructions
Zc <- function(AAcomp, ...) {
  # The number of carbon atoms in each amino acid
  nC_AA <- c(Ala = 3, Cys = 3, Asp = 4, Glu = 5, Phe = 9, Gly = 2, His = 6, 
    Ile = 6, Lys = 6, Leu = 6, Met = 5, Asn = 4, Pro = 5, Gln = 5, 
    Arg = 6, Ser = 3, Thr = 4, Val = 5, Trp = 11, Tyr = 9)
  # The Ztot of the amino acids == CHNOSZ::ZC(info(info(aminoacids("")))$formula) * nC_AA
  Ztot_AA <- c(Ala = 0, Cys = 2, Asp = 4, Glu = 2, Phe = -4, Gly = 2, His = 4, 
    Ile = -6, Lys = -4, Leu = -6, Met = -2, Asn = 4, Pro = -2, Gln = 2, 
    Arg = 2, Ser = 2, Thr = 0, Val = -4, Trp = -2, Tyr = -2)
  # The Zc of the amino acids == CHNOSZ::ZC(info(info(aminoacids("")))$formula)
  Zc_AA <- Ztot_AA / nC_AA
  # Find columns with names for the amino acids
  isAA <- tolower(colnames(AAcomp)) %in% tolower(names(Zc_AA))
  iAA <- match(tolower(colnames(AAcomp)[isAA]), tolower(names(Zc_AA)))
  # Calculate the nC for all occurrences of each amino acid
  multC <- t(t(AAcomp[, isAA, drop = FALSE]) * nC_AA[iAA])
  # Multiply nC by Zc
  multZc <- t(t(multC) * Zc_AA[iAA])
  # Calculate the total Zc and nC, then the overall Zc
  Zctot <- rowSums(multZc)
  nCtot <- rowSums(multC)
  Zctot / nCtot
}

# Per-residue stoichiometric hydration state 20181228
# Add 'terminal_H2O' argument 20221018
nH2O <- function(AAcomp, basis = "QEC", terminal_H2O = 0) {
  if(basis == "QEC") {
    # To get the number of H2O in reactions to form amino acid residues from the "QEC" basis:
    ## library(CHNOSZ)
    ## basis("QEC")
    ## nH2O_AA <- species(aminoacids(""))$H2O
    ## names(nH2O_AA) <- aminoacids(3)
    nH2O_AA <- c( Ala =  0.6, Cys =    0, Asp = -0.2, Glu =    0, Phe = -2.2, Gly =  0.4, His = -1.8,
      Ile =  1.2, Lys =  1.2, Leu =  1.2, Met =  0.4, Asn = -0.2, Pro =    0, Gln =    0,
      Arg =  0.2, Ser =  0.6, Thr =  0.8, Val =    1, Trp = -3.8, Tyr = -2.2) - 1
    # Note: subtraction of 1 is to get amino acid residues
  }
  # QCa basis species 20200818
  if(basis == "QCa") {
    ## library(CHNOSZ)
    ## basis(c("cysteine", "glutamine", "acetic acid", "H2O", "O2"))
    ## nH2O_AA <- species(aminoacids(""))$H2O
    ## names(nH2O_AA) <- aminoacids(3)
    nH2O_AA <- c(Ala = 0.5, Cys = 0, Asp = -0.5, Glu = -0.5, Phe = -3.5, Gly = 0.5,
      His = -1.5, Ile = 0.5, Lys = 1, Leu = 0.5, Met = 0, Asn = 0,
      Pro = -0.5, Gln = 0, Arg = 1, Ser = 0.5, Thr = 0.5, Val = 0.5,
      Trp = -5, Tyr = -3.5) - 1
  }
  # Find columns with names for the amino acids
  isAA <- tolower(colnames(AAcomp)) %in% tolower(names(nH2O_AA))
  iAA <- match(tolower(colnames(AAcomp)[isAA]), tolower(names(nH2O_AA)))
  # Calculate total number of H2O in reactions to form proteins
  nH2O <- rowSums(t(t(AAcomp[, isAA, drop = FALSE]) * nH2O_AA[iAA]))
  # Add one H2O for each pair of terminal groups (i.e., number of polypeptide chains)
  nH2O <- nH2O + terminal_H2O
  # Divide by number of residues (length of protein)
  nH2O / rowSums(AAcomp[, isAA, drop = FALSE])
}

# Per-residue stoichiometric oxidation state 20201016
nO2 <- function(AAcomp, basis = "QEC", ...) {
  if(basis == "QEC") {
    # To get the number of O2 in reactions to form amino acid residues from the "QEC" basis:
    ## library(CHNOSZ)
    ## basis("QEC")
    ## nO2_AA <- species(aminoacids(""))[["O2"]]
    ## names(nO2_AA) <- aminoacids(3)
    nO2_AA <- c(Ala = -0.3, Cys = 0, Asp = 0.6, Glu = 0, Phe = -1.9, Gly = 0.3, 
      His = 0.4, Ile = -2.1, Lys = -1.6, Leu = -2.1, Met = -1.2, Asn = 0.6, 
      Pro = -1, Gln = 0, Arg = -0.1, Ser = 0.2, Thr = -0.4, Val = -1.5, 
      Trp = -1.6, Tyr = -1.4)
  }
  # QCa basis species 20200818
  if(basis == "QCa") {
    ## library(CHNOSZ)
    ## basis(c("cysteine", "glutamine", "acetic acid", "H2O", "O2"))
    ## nO2_AA <- species(aminoacids(""))$O2
    ## names(nO2_AA) <- aminoacids(3)
    nO2_AA <- c(Ala = -0.25, Cys = 0, Asp = 0.75, Glu = 0.25, Phe = -1.25, 
      Gly = 0.25, His = 0.25, Ile = -1.75, Lys = -1.5, Leu = -1.75, 
      Met = -1, Asn = 0.5, Pro = -0.75, Gln = 0, Arg = -0.5, Ser = 0.25, 
      Thr = -0.25, Val = -1.25, Trp = -1, Tyr = -0.75)
  }
  # Find columns with names for the amino acids
  isAA <- tolower(colnames(AAcomp)) %in% tolower(names(nO2_AA))
  iAA <- match(tolower(colnames(AAcomp)[isAA]), tolower(names(nO2_AA)))
  # Calculate total number of O2 in reactions to form proteins
  nO2 <- rowSums(t(t(AAcomp[, isAA, drop = FALSE]) * nO2_AA[iAA]))
  # Divide by number of residues (length of protein)
  nO2 / rowSums(AAcomp[, isAA, drop = FALSE])
}

# Grand average of hydropathy (GRAVY) 20191024
GRAVY <- function(AAcomp, ...) {
  # Values of the hydropathy index from Kyte and Doolittle, 1982
  # doi:10.1016/0022-2836(82)90515-0
  Hind <- c(Ala =  1.8, Cys =  2.5, Asp = -3.5, Glu = -3.5, Phe =  2.8,
            Gly = -0.4, His = -3.2, Ile =  4.5, Lys = -3.9, Leu =  3.8,
            Met =  1.9, Asn = -3.5, Pro = -1.6, Gln = -3.5, Arg = -4.5,
            Ser = -0.8, Thr = -0.7, Val =  4.2, Trp = -0.9, Tyr = -1.3)
  # Find columns with names for the amino acids
  isAA <- colnames(AAcomp) %in% names(Hind)
  iAA <- match(colnames(AAcomp)[isAA], names(Hind))
  # Calculate total of hydropathy values for each protein
  sumHind <- rowSums(t(t(AAcomp[, isAA, drop = FALSE]) * Hind[iAA]))
  # Divide by length of proteins to get GRAVY
  sumHind / rowSums(AAcomp[, isAA, drop = FALSE])
}

# Isoelectric point 20191026
pI <- function(AAcomp, terminal_H2O = 1, ...) {
  if(nrow(AAcomp) == 0) return(numeric(0))
  # A function to calculate isoelectric point for a single amino acid composition
  onepI <- function(AA) {
    # Find the column names of AAcomp that are in Ztab
    isZ <- names(AA) %in% dimnames(Ztab)[[2]]
    iZ <- match(names(AA)[isZ], dimnames(Ztab)[[2]])
    # Calculate the total charge as a function of pH
    # ... the "else" is in case we have a data frame (used when first writing this function)
    if(is.numeric(AA)) Ztot <- Ztab[, iZ] %*% AA[isZ] else
    Ztot <- Ztab[, iZ] %*% as.matrix(t(AA[, isZ]))
    # Find pH where charge is closest to zero
    # (absolute charge is minimized)
    ipH <- which.min(abs(Ztot))
    pI <- Ztab[ipH, 1]
    if(length(pI) == 0) pI <- NA
    pI
  }
  # Number of N- and C-terminal groups
  Nterm <- Cterm <- terminal_H2O
  if(!"Nterm" %in% names(AAcomp)) AAcomp <- cbind(AAcomp, Nterm = Nterm)
  if(!"Cterm" %in% names(AAcomp)) AAcomp <- cbind(AAcomp, Cterm = Cterm)
  # NOTE: apply() converts the input to matrix,
  # so we extract the numeric columns of AAcomp to avoid possible coercion of all values to character
  isnum <- unlist(lapply(AAcomp, "class")) %in% c("integer", "numeric")
  myAA <- AAcomp[, isnum, drop = FALSE]
  # Run the calculation for each composition
  apply(myAA, 1, onepI)
}

# Per-residue molecular weight 20200501
MW <- function(AAcomp, terminal_H2O = 0, ...) {
  # Mass per residue:
  # MW_AA <- sapply(CHNOSZ::makeup(info(aminoacids(""))), mass) - mass("H2O")
  # names(MW_AA) <- aminoacids(3)
  MW_AA <- c(Ala = 71.0788, Cys = 103.1388, Asp = 115.0886, Glu = 129.11548, 
    Phe = 147.17656, Gly = 57.05192, His = 137.14108, Ile = 113.15944, 
    Lys = 128.17408, Leu = 113.15944, Met = 131.19256, Asn = 114.10384, 
    Pro = 97.11668, Gln = 128.13072, Arg = 156.18748, Ser = 87.0782, 
    Thr = 101.10508, Val = 99.13256, Trp = 186.2132, Tyr = 163.17596
  )
  # Find columns with names for the amino acids
  isAA <- tolower(colnames(AAcomp)) %in% tolower(names(MW_AA))
  iAA <- match(colnames(AAcomp)[isAA], names(MW_AA))
  # Calculate total MW of residues in each protein
  MW <- rowSums(t(t(AAcomp[, isAA, drop = FALSE]) * MW_AA[iAA]))
  # Add mass of H2O for each pair of terminal groups
  # MW_H2O <- mass("H2O")
  MW_H2O <- 18.01528
  MW <- MW + terminal_H2O * MW_H2O
  # Divide by number of residues (length of protein)
  MW / rowSums(AAcomp[, isAA, drop = FALSE])
}

# Per-protein molecular weight 20240304
pMW <- function(AAcomp, terminal_H2O = 1, ...) {
  # Mass per residue:
  # MW_AA <- sapply(CHNOSZ::makeup(info(aminoacids(""))), mass) - mass("H2O")
  # names(MW_AA) <- aminoacids(3)
  MW_AA <- c(Ala = 71.0788, Cys = 103.1388, Asp = 115.0886, Glu = 129.11548, 
    Phe = 147.17656, Gly = 57.05192, His = 137.14108, Ile = 113.15944, 
    Lys = 128.17408, Leu = 113.15944, Met = 131.19256, Asn = 114.10384, 
    Pro = 97.11668, Gln = 128.13072, Arg = 156.18748, Ser = 87.0782, 
    Thr = 101.10508, Val = 99.13256, Trp = 186.2132, Tyr = 163.17596
  )
  # Find columns with names for the amino acids
  isAA <- tolower(colnames(AAcomp)) %in% tolower(names(MW_AA))
  iAA <- match(colnames(AAcomp)[isAA], names(MW_AA))
  # Calculate total MW of residues in each protein
  MW <- rowSums(t(t(AAcomp[, isAA, drop = FALSE]) * MW_AA[iAA]))
  # Add mass of H2O for each pair of terminal groups
  # MW_H2O <- mass("H2O")
  MW_H2O <- 18.01528
  MW + terminal_H2O * MW_H2O
}

# Per-residue volume 20240301
V0 <- function(AAcomp, terminal_H2O = 0, ...) {
  # Volume per residue using group contributions from Dick et al., 2006:
  # i.e. residue = [sidechain group] + [backbone group]
  # V0_AA <- info(info(paste0("[", aminoacids(3), "]")))$V + info(info("[UPBB]"))$V
  # names(V0_AA) <- aminoacids(3)
  V0_AA <- c(Ala = 53.16, Cys = 66.209, Asp = 67.412, Glu = 82.917, Phe = 114.841, 
    Gly = 35.902, His = 92.049, Ile = 98.5, Lys = 101.344, Leu = 100.496, 
    Met = 98.128, Asn = 70.122, Pro = 75.345, Gln = 86.374, Arg = 132.121, 
    Ser = 53.338, Thr = 70.326, Val = 83.575, Trp = 136.341, Tyr = 117.2
  )
  # Find columns with names for the amino acids
  isAA <- tolower(colnames(AAcomp)) %in% tolower(names(V0_AA))
  iAA <- match(colnames(AAcomp)[isAA], names(V0_AA))
  # Calculate total V0 of residues in each protein
  V0 <- rowSums(t(t(AAcomp[, isAA, drop = FALSE]) * V0_AA[iAA]))
  # Add volume of H2O for each pair of terminal groups
  # V0_H2O <- info(info("[AABB]"))$V - info(info("[UPBB]"))$V
  V0_H2O <- 7.289
  V0 <- V0 + terminal_H2O * V0_H2O
  # Divide by number of residues (length of protein)
  V0 / rowSums(AAcomp[, isAA, drop = FALSE])
}

# Per-protein volume 20240301
pV0 <- function(AAcomp, terminal_H2O = 1, ...) {
  # Volume per residue using group contributions from Dick et al., 2006:
  # i.e. residue = [sidechain group] + [backbone group]
  # V0_AA <- info(info(paste0("[", aminoacids(3), "]")))$V + info(info("[UPBB]"))$V
  # names(V0_AA) <- aminoacids(3)
  V0_AA <- c(Ala = 53.16, Cys = 66.209, Asp = 67.412, Glu = 82.917, Phe = 114.841, 
    Gly = 35.902, His = 92.049, Ile = 98.5, Lys = 101.344, Leu = 100.496, 
    Met = 98.128, Asn = 70.122, Pro = 75.345, Gln = 86.374, Arg = 132.121, 
    Ser = 53.338, Thr = 70.326, Val = 83.575, Trp = 136.341, Tyr = 117.2
  )
  # Find columns with names for the amino acids
  isAA <- tolower(colnames(AAcomp)) %in% tolower(names(V0_AA))
  iAA <- match(colnames(AAcomp)[isAA], names(V0_AA))
  # Calculate total V0 of residues in each protein
  V0 <- rowSums(t(t(AAcomp[, isAA, drop = FALSE]) * V0_AA[iAA]))
  # Add volume of H2O for each pair of terminal groups
  # V0_H2O <- info(info("[AABB]"))$V - info(info("[UPBB]"))$V
  V0_H2O <- 7.289
  V0 + terminal_H2O * V0_H2O
}

# Per-residue entropy 20240308
S0 <- function(AAcomp, terminal_H2O = 0, ...) {
  # Entropy per residue using group contributions from Dick et al., 2006:
  # i.e. residue = [sidechain group] + [backbone group]
  # S0_AA <- info(info(paste0("[", aminoacids(3), "]")))$S + info(info("[UPBB]"))$S
  # names(S0_AA) <- aminoacids(3)
  S0_AA <- c(Ala = 19.86, Cys = 27.35, Asp = 36.25, Glu = 42.23, Phe = 37.63, 
    Gly = 18.92, His = 47.03, Ile = 30.73, Lys = 42.31, Leu = 31.44, 
    Met = 43.83, Asn = 38.91, Pro = 30.86, Gln = 43.44, Arg = 57.76, 
    Ser = 28.27, Thr = 25.26, Val = 26.71, Trp = 40.99, Tyr = 40.44
  )
  # Find columns with names for the amino acids
  isAA <- tolower(colnames(AAcomp)) %in% tolower(names(S0_AA))
  iAA <- match(colnames(AAcomp)[isAA], names(S0_AA))
  # Calculate total S0 of residues in each protein
  S0 <- rowSums(t(t(AAcomp[, isAA, drop = FALSE]) * S0_AA[iAA]))
  # Add volume of H2O for each pair of terminal groups
  # S0_H2O <- info(info("[AABB]"))$S - info(info("[UPBB]"))$S
  S0_H2O <- 18.97
  S0 <- S0 + terminal_H2O * S0_H2O
  # Divide by number of residues (length of protein)
  S0 / rowSums(AAcomp[, isAA, drop = FALSE])
}

# Per-protein entropy 20240308
pS0 <- function(AAcomp, terminal_H2O = 1, ...) {
  # Entropy per residue using group contributions from Dick et al., 2006:
  # i.e. residue = [sidechain group] + [backbone group]
  # S0_AA <- info(info(paste0("[", aminoacids(3), "]")))$S + info(info("[UPBB]"))$S
  # names(S0_AA) <- aminoacids(3)
  S0_AA <- c(Ala = 19.86, Cys = 27.35, Asp = 36.25, Glu = 42.23, Phe = 37.63, 
    Gly = 18.92, His = 47.03, Ile = 30.73, Lys = 42.31, Leu = 31.44, 
    Met = 43.83, Asn = 38.91, Pro = 30.86, Gln = 43.44, Arg = 57.76, 
    Ser = 28.27, Thr = 25.26, Val = 26.71, Trp = 40.99, Tyr = 40.44
  )
  # Find columns with names for the amino acids
  isAA <- tolower(colnames(AAcomp)) %in% tolower(names(S0_AA))
  iAA <- match(colnames(AAcomp)[isAA], names(S0_AA))
  # Calculate total S0 of residues in each protein
  S0 <- rowSums(t(t(AAcomp[, isAA, drop = FALSE]) * S0_AA[iAA]))
  # Add volume of H2O for each pair of terminal groups
  # S0_H2O <- info(info("[AABB]"))$S - info(info("[UPBB]"))$S
  S0_H2O <- 18.97
  S0 + terminal_H2O * S0_H2O
}

# Protein length 20200501
plength <- function(AAcomp, ...) {
  AA_names <- c(
    "Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu", "Met",
    "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr"
  )
  # Find columns with names for the amino acids
  isAA <- tolower(colnames(AAcomp)) %in% tolower(AA_names)
  # Sum amino acid counts to get protein length
  rowSums(AAcomp[, isAA])
}

# H/C (H:C ratio) 20230707
HC <- function(AAcomp, ...) {
  # The number of H in each amino acid residue; calculated in CHNOSZ:
  # nH_AA <- sapply(makeup(info(info(aminoacids("")))$formula), "[", "H")
  # nH_AA <- nH_AA - 2  # Take H-OH off of amino acids to make residues
  # names(nH_AA) <- aminoacids(3)
  nH_AA <- c(Ala = 5, Cys = 5, Asp = 5, Glu = 7, Phe = 9, Gly = 3, His = 7,
    Ile = 11, Lys = 12, Leu = 11, Met = 9, Asn = 6, Pro = 7, Gln = 8,
    Arg = 12, Ser = 5, Thr = 7, Val = 9, Trp = 10, Tyr = 9)
  # The number of carbon atoms in each amino acid
  nC_AA <- c(Ala = 3, Cys = 3, Asp = 4, Glu = 5, Phe = 9, Gly = 2, His = 6, 
    Ile = 6, Lys = 6, Leu = 6, Met = 5, Asn = 4, Pro = 5, Gln = 5, 
    Arg = 6, Ser = 3, Thr = 4, Val = 5, Trp = 11, Tyr = 9)
  # Find columns with names for the amino acids
  isAA <- tolower(colnames(AAcomp)) %in% tolower(names(nH_AA))
  iAA <- match(tolower(colnames(AAcomp)[isAA]), tolower(names(nH_AA)))
  # Count the number of C in all residues
  numC <- t(t(AAcomp[, isAA, drop = FALSE]) * nC_AA[iAA])
  # Count the number of H in all residues
  numH <- t(t(AAcomp[, isAA, drop = FALSE]) * nH_AA[iAA])
  # Calculate the total number of H and C, then the overall H/C
  Htot <- rowSums(numH)
  Ctot <- rowSums(numC)
  Htot / Ctot
}

# N/C (N:C ratio) 20230707
NC <- function(AAcomp, ...) {
  # The number of N in each amino acid residue; calculated in CHNOSZ:
  # nN_AA <- sapply(makeup(info(info(aminoacids("")))$formula), "[", "N")
  # names(nN_AA) <- aminoacids(3)
  nN_AA <- c(Ala = 1, Cys = 1, Asp = 1, Glu = 1, Phe = 1, Gly = 1, His = 3,
  Ile = 1, Lys = 2, Leu = 1, Met = 1, Asn = 2, Pro = 1, Gln = 2,
  Arg = 4, Ser = 1, Thr = 1, Val = 1, Trp = 2, Tyr = 1)
  # The number of carbon atoms in each amino acid
  nC_AA <- c(Ala = 3, Cys = 3, Asp = 4, Glu = 5, Phe = 9, Gly = 2, His = 6, 
    Ile = 6, Lys = 6, Leu = 6, Met = 5, Asn = 4, Pro = 5, Gln = 5, 
    Arg = 6, Ser = 3, Thr = 4, Val = 5, Trp = 11, Tyr = 9)
  # Find columns with names for the amino acids
  isAA <- tolower(colnames(AAcomp)) %in% tolower(names(nN_AA))
  iAA <- match(tolower(colnames(AAcomp)[isAA]), tolower(names(nN_AA)))
  # Count the number of C in all residues
  numC <- t(t(AAcomp[, isAA, drop = FALSE]) * nC_AA[iAA])
  # Count the number of N in all residues
  numN <- t(t(AAcomp[, isAA, drop = FALSE]) * nN_AA[iAA])
  # Calculate the total number of N and C, then the overall N/C
  Ntot <- rowSums(numN)
  Ctot <- rowSums(numC)
  Ntot / Ctot
}

# O/C (O:C ratio) 20230707
OC <- function(AAcomp, ...) {
  # The number of O in each amino acid residue; calculated in CHNOSZ:
  # nO_AA <- sapply(makeup(info(info(aminoacids("")))$formula), "[", "O")
  # nO_AA <- nO_AA - 1  # Take H-OH off of amino acids to make residues
  # names(nO_AA) <- aminoacids(3)
  nO_AA <- c(Ala = 1, Cys = 1, Asp = 3, Glu = 3, Phe = 1, Gly = 1, His = 1,
  Ile = 1, Lys = 1, Leu = 1, Met = 1, Asn = 2, Pro = 1, Gln = 2,
  Arg = 1, Ser = 2, Thr = 2, Val = 1, Trp = 1, Tyr = 2)
  # The number of carbon atoms in each amino acid
  nC_AA <- c(Ala = 3, Cys = 3, Asp = 4, Glu = 5, Phe = 9, Gly = 2, His = 6, 
    Ile = 6, Lys = 6, Leu = 6, Met = 5, Asn = 4, Pro = 5, Gln = 5, 
    Arg = 6, Ser = 3, Thr = 4, Val = 5, Trp = 11, Tyr = 9)
  # Find columns with names for the amino acids
  isAA <- tolower(colnames(AAcomp)) %in% tolower(names(nO_AA))
  iAA <- match(tolower(colnames(AAcomp)[isAA]), tolower(names(nO_AA)))
  # Count the number of C in all residues
  numC <- t(t(AAcomp[, isAA, drop = FALSE]) * nC_AA[iAA])
  # Count the number of O in all residues
  numO <- t(t(AAcomp[, isAA, drop = FALSE]) * nO_AA[iAA])
  # Calculate the total number of O and C, then the overall O/C
  Otot <- rowSums(numO)
  Ctot <- rowSums(numC)
  Otot / Ctot
}

# S/C (S:C ratio) 20230707
SC <- function(AAcomp, ...) {
  # The number of S in each amino acid residue; calculated in CHNOSZ:
  # nS_AA <- sapply(makeup(info(info(aminoacids("")))$formula), "[", "S")
  # nS_AA[is.na(nS_AA)] <- 0
  # names(nS_AA) <- aminoacids(3)
  nS_AA <- c(Ala = 0, Cys = 1, Asp = 0, Glu = 0, Phe = 0, Gly = 0, His = 0,
  Ile = 0, Lys = 0, Leu = 0, Met = 1, Asn = 0, Pro = 0, Gln = 0,
  Arg = 0, Ser = 0, Thr = 0, Val = 0, Trp = 0, Tyr = 0)
  # The number of carbon atoms in each amino acid
  nC_AA <- c(Ala = 3, Cys = 3, Asp = 4, Glu = 5, Phe = 9, Gly = 2, His = 6, 
    Ile = 6, Lys = 6, Leu = 6, Met = 5, Asn = 4, Pro = 5, Gln = 5, 
    Arg = 6, Ser = 3, Thr = 4, Val = 5, Trp = 11, Tyr = 9)
  # Find columns with names for the amino acids
  isAA <- tolower(colnames(AAcomp)) %in% tolower(names(nS_AA))
  iAA <- match(tolower(colnames(AAcomp)[isAA]), tolower(names(nS_AA)))
  # Count the number of C in all residues
  numC <- t(t(AAcomp[, isAA, drop = FALSE]) * nC_AA[iAA])
  # Count the number of S in all residues
  numS <- t(t(AAcomp[, isAA, drop = FALSE]) * nS_AA[iAA])
  # Calculate the total number of S and C, then the overall S/C
  Stot <- rowSums(numS)
  Ctot <- rowSums(numC)
  Stot / Ctot
}

# Number of carbon atoms per residue 20240322
nC <- function(AAcomp, ...) {
  # The number of carbon atoms in each amino acid
  nC_AA <- c(Ala = 3, Cys = 3, Asp = 4, Glu = 5, Phe = 9, Gly = 2, His = 6, 
    Ile = 6, Lys = 6, Leu = 6, Met = 5, Asn = 4, Pro = 5, Gln = 5, 
    Arg = 6, Ser = 3, Thr = 4, Val = 5, Trp = 11, Tyr = 9)
  # Find columns with names for the amino acids
  isAA <- tolower(colnames(AAcomp)) %in% tolower(names(nC_AA))
  iAA <- match(tolower(colnames(AAcomp)[isAA]), tolower(names(nC_AA)))
  # Calculate the nC for all occurrences of each amino acid
  multC <- t(t(AAcomp[, isAA, drop = FALSE]) * nC_AA[iAA])
  # Calculate the total nC and divide by number of residues (length of protein)
  rowSums(multC) / rowSums(AAcomp[, isAA, drop = FALSE])
}

# Number of carbon atoms per protein 20240322
pnC <- function(AAcomp, ...) {
  # The number of carbon atoms in each amino acid
  nC_AA <- c(Ala = 3, Cys = 3, Asp = 4, Glu = 5, Phe = 9, Gly = 2, His = 6, 
    Ile = 6, Lys = 6, Leu = 6, Met = 5, Asn = 4, Pro = 5, Gln = 5, 
    Arg = 6, Ser = 3, Thr = 4, Val = 5, Trp = 11, Tyr = 9)
  # Find columns with names for the amino acids
  isAA <- tolower(colnames(AAcomp)) %in% tolower(names(nC_AA))
  iAA <- match(tolower(colnames(AAcomp)[isAA]), tolower(names(nC_AA)))
  # Calculate the nC for all occurrences of each amino acid
  multC <- t(t(AAcomp[, isAA, drop = FALSE]) * nC_AA[iAA])
  # Calculate the total nC
  rowSums(multC)
}

# Density 20240304
Density <- function(AAcomp, ...) {
  MW(AAcomp, ...) / V0(AAcomp, ...)
}

# Specific volume 20240308
V0g <- function(AAcomp, ...) {
  V0(AAcomp, ...) / MW(AAcomp, ...)
}

# Specific entropy 20240308
S0g <- function(AAcomp, ...) {
  S0(AAcomp, ...) / MW(AAcomp, ...)
}

# Entropy density 20240309
SV <- function(AAcomp, ...) {
  S0(AAcomp, ...) / V0(AAcomp, ...)
}

# Specific Zc 20240317
Zcg <- function(AAcomp, ...) {
  Zc(AAcomp, ...) * nC(AAcomp, ...) / MW(AAcomp, ...)
}

# Specific nO2 20240317
nO2g <- function(AAcomp, ...) {
  nO2(AAcomp, ...) / MW(AAcomp, ...)
}

# Specific nH2O 20240317
nH2Og <- function(AAcomp, ...) {
  nH2O(AAcomp, ...) / MW(AAcomp, ...)
}


# Metabolic cost 20211220
# Moved from JMDplots/evdevH2O.R to canprot 20240304
Cost <- function(AAcomp, ...) {
  # Amino acid cost from Akashi and Gojobori (2002)
  # doi:10.1073/pnas.062526999
  Cost_AA <- c(Ala = 11.7, Cys = 24.7, Asp = 12.7, Glu = 15.3, Phe = 52.0,
               Gly = 11.7, His = 38.3, Ile = 32.3, Lys = 30.3, Leu = 27.3,
               Met = 34.3, Asn = 14.7, Pro = 20.3, Gln = 16.3, Arg = 27.3,
               Ser = 11.7, Thr = 18.7, Val = 23.3, Trp = 74.3, Tyr = 50.0)
  # Find columns with names for the amino acids
  isAA <- tolower(colnames(AAcomp)) %in% tolower(names(Cost_AA))
  iAA <- match(colnames(AAcomp)[isAA], names(Cost_AA))
  # Calculate total Cost of residues in each protein
  Cost <- rowSums(t(t(AAcomp[, isAA, drop = FALSE]) * Cost_AA[iAA]))
  # Divide by number of residues (length of protein)
  Cost / rowSums(AAcomp[, isAA, drop = FALSE])
}

# Respiratory cost 20240304
RespiratoryCost <- function(AAcomp, ...) {
  # Respiratory cost from Wagner (2005)
  # doi:10.1093/molbev/msi126
  Cost_AA <- c(Ala = 14.5, Cys = 26.5, Asp = 15.5, Glu =  9.5, Phe = 61.0,
               Gly = 14.5, His = 29.0, Ile = 38.0, Lys = 36.0, Leu = 37.0,
               Met = 36.5, Asn = 18.5, Pro = 14.5, Gln = 10.5, Arg = 20.5,
               Ser = 14.5, Thr = 21.5, Val = 29.0, Trp = 75.5, Tyr = 59.0)
  # Find columns with names for the amino acids
  isAA <- tolower(colnames(AAcomp)) %in% tolower(names(Cost_AA))
  iAA <- match(colnames(AAcomp)[isAA], names(Cost_AA))
  # Calculate total Cost of residues in each protein
  Cost <- rowSums(t(t(AAcomp[, isAA, drop = FALSE]) * Cost_AA[iAA]))
  # Divide by number of residues (length of protein)
  Cost / rowSums(AAcomp[, isAA, drop = FALSE])
}

# Fermentative cost 20240304
FermentativeCost <- function(AAcomp, ...) {
  # Fermentative cost from Wagner (2005)
  # doi:10.1093/molbev/msi126
  Cost_AA <- c(Ala =  2, Cys = 13, Asp =  3, Glu =  2, Phe = 10,
               Gly =  1, His =  5, Ile = 14, Lys = 12, Leu =  4,
               Met = 24, Asn =  6, Pro =  7, Gln =  3, Arg = 13,
               Ser =  1, Thr =  9, Val =  4, Trp = 14, Tyr =  8)
  # Find columns with names for the amino acids
  isAA <- tolower(colnames(AAcomp)) %in% tolower(names(Cost_AA))
  iAA <- match(colnames(AAcomp)[isAA], names(Cost_AA))
  # Calculate total Cost of residues in each protein
  Cost <- rowSums(t(t(AAcomp[, isAA, drop = FALSE]) * Cost_AA[iAA]))
  # Divide by number of residues (length of protein)
  Cost / rowSums(AAcomp[, isAA, drop = FALSE])
}

# Biosynthetic cost in bacteria 20240304
B20Cost <- function(AAcomp, ...) {
  # Biosynthetic cost from Zhang et al. (2018)
  # doi:10.1038/s41467-018-06461-1
  Cost_AA <- c(Ala =  11.7, Cys = 741.0, Asp = 114.3, Glu =  76.5, Phe = 208.0,
               Gly =  11.7, His = 536.2, Ile =  64.6, Lys = 242.4, Leu =  54.6,
               Met = 445.9, Asn = 147.0, Pro =  60.9, Gln = 130.4, Arg = 109.2,
               Ser =  70.2, Thr = 112.2, Val =  46.6, Trp = 891.6, Tyr = 350.0)
  # Find columns with names for the amino acids
  isAA <- tolower(colnames(AAcomp)) %in% tolower(names(Cost_AA))
  iAA <- match(colnames(AAcomp)[isAA], names(Cost_AA))
  # Calculate total Cost of residues in each protein
  Cost <- rowSums(t(t(AAcomp[, isAA, drop = FALSE]) * Cost_AA[iAA]))
  # Divide by number of residues (length of protein)
  Cost / rowSums(AAcomp[, isAA, drop = FALSE])
}

# Biosynthetic cost in yeast 20240304
Y20Cost <- function(AAcomp, ...) {
  # Biosynthetic cost from Zhang et al. (2018)
  # doi:10.1038/s41467-018-06461-1
  Cost_AA <- c(Ala =  14.5, Cys = 795.0, Asp = 139.5, Glu =  47.5, Phe = 244.0,
               Gly =  14.5, His = 406.0, Ile =  76.0, Lys = 288.0, Leu =  74.0,
               Met = 474.5, Asn = 185.0, Pro =  43.5, Gln =  84.0, Arg =  82.0,
               Ser =  87.0, Thr = 129.0, Val =  58.0, Trp = 906.0, Tyr = 413.0)
  # Find columns with names for the amino acids
  isAA <- tolower(colnames(AAcomp)) %in% tolower(names(Cost_AA))
  iAA <- match(colnames(AAcomp)[isAA], names(Cost_AA))
  # Calculate total Cost of residues in each protein
  Cost <- rowSums(t(t(AAcomp[, isAA, drop = FALSE]) * Cost_AA[iAA]))
  # Divide by number of residues (length of protein)
  Cost / rowSums(AAcomp[, isAA, drop = FALSE])
}

# Biosynthetic cost in humans 20240304
H11Cost <- function(AAcomp, ...) {
  # Biosynthetic cost from Zhang et al. (2018)
  # doi:10.1038/s41467-018-06461-1
  Cost_AA <- c(Ala =  12.5, Cys = 330.0, Asp = 121.5, Glu =  42.5, Phe =   0.0,
               Gly =  11.0, His =   0.0, Ile =   0.0, Lys =   0.0, Leu =   0.0,
               Met =   0.0, Asn = 165.0, Pro =  43.5, Gln =  76.0, Arg =  58.0,
               Ser =  66.0, Thr =   0.0, Val =   0.0, Trp =   0.0, Tyr =  17.5)
  # Find columns with names for the amino acids
  isAA <- tolower(colnames(AAcomp)) %in% tolower(names(Cost_AA))
  iAA <- match(colnames(AAcomp)[isAA], names(Cost_AA))
  # Calculate total Cost of residues in each protein
  Cost <- rowSums(t(t(AAcomp[, isAA, drop = FALSE]) * Cost_AA[iAA]))
  # Divide by number of residues (length of protein)
  Cost / rowSums(AAcomp[, isAA, drop = FALSE])
}

# Now that we've added the functions, check that they are all listed in cplab
cplist1 <- ls()
metrics <- setdiff(cplist1, c(cplist0, "cplist0"))
not_in_cplab <- metrics[!metrics %in% names(cplab)]
if(length(not_in_cplab) > 0) stop(paste("These metrics are not in cplab:", paste(not_in_cplab, collapse = ", ")))
not_in_metrics <- names(cplab)[!names(cplab) %in% metrics]
if(length(not_in_metrics) > 0) stop(paste("These metrics in cplab don't have functions:", paste(not_in_metrics, collapse = ", ")))

#########################
### UNEXPORTED OBJECT ###
### ( used in pI() )  ###
#########################

# Tabulate charges for sidechains and terminal groups from pH 0 to 14
Ztab <- local({
  # A function to calculate charge as a function of pH for a single group
  ZpH <- function(pK, Z, pH) {
    alpha <- 1/(1 + 10^(Z * (pH - pK)))
    alpha * Z
  }
  # List the pKs of the groups
  pK <- list(Cterm = 3.55, Nterm = 7.5,
    Asp = 4.05, Glu = 4.45, His = 5.98,
    Cys = 9, Lys = 10, Tyr = 10, Arg = 12
  )
  # List the unit charges of the groups
  Z <- list(Cterm = -1, Nterm = 1,
    Asp = -1, Glu = -1, His = 1,
    Cys = -1, Lys = 1, Tyr = -1, Arg = 1
  )
  # Get the charges for a range of pH values
  pH <- seq(0, 14, 0.01)
  Ztab <- mapply(ZpH, pK = pK, Z = Z, MoreArgs = list(pH = pH))
  # Add a column with the pH values
  cbind(pH = pH, Ztab)
})
