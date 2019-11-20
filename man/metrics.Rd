\encoding{UTF-8}
\name{metrics}
\alias{metrics}
\alias{H2OAA}
\alias{ZCAA}
\alias{GRAVY}
\alias{pI}
\title{Calculate Compositional Metrics for Proteins}
\description{
These functions calculate compositional metrics of proteins given a data frame of amino acid compositions.
}

\usage{
  ZCAA(AAcomp, nothing = NULL)
  H2OAA(AAcomp, basis = "rQEC")
  GRAVY(AAcomp)
  pI(AAcomp)
}

\arguments{
  \item{AAcomp}{data frame, amino acid compositions}
  \item{nothing}{dummy argument}
  \item{basis}{character, basis species}
}

\details{
\code{ZCAA} calculates average oxidation state of carbon (\ZC).
\code{H2OAA} calculates stoichiometric hydration state (\nH2O) per residue.
\code{GRAVY} calculates the grand average of hydropathicity.
\code{pI} calculates the isoelectric point.

These functions calculate the metrics from a data frame of amino acid compositions.
Columns in \code{AAcomp} should be named with the three-letter abbreviations for the amino acids (e.g. \samp{Ala}).
For calculations of \ZC from chemical formulas of organic molecules (not only proteins), use the \code{\link[CHNOSZ]{ZC}} function in \pkg{CHNOSZ}.

The default for \code{basis} stands for relative water content, calculated as residuals of a linear model for \nH2O vs \ZC of the 20 common amino acids using the basis species glutamine - glutamic acid - cysteine - \H2O - \O2.
Furthermore, a constant of 0.355 is subtracted from these residuals to make the mean for all human proteins = 0.
This can be changed to \samp{QEC} to use the stoichiometric water content calculated directly from the basis species (as in Dick, 2017).
See Dick et al. (2019) for more information.

Note that \ZC is a per-carbon average, but \nH2O is a per-residue average.
The contribution of \H2O from the terminal groups of proteins is counted, so shorter proteins have slightly greater \nH2O.

\code{nothing} is an extra argument that does nothing.
It is provided so that these functions can be called with the same numbers of arguments.

For \code{GRAVY}, values of the hydropathy index for individual amino acids are from Kyte and Doolittle (1982).

For \code{pI}, the net charge for each ionizable group is pre-calculated from pH 0 to 14 at intervals of 0.01.
The isoelectric point is found as the pH where the sum of charges of all groups in the protein is closest to zero.
The \pK values for the terminal groups and sidechains are taken from Bjellqvist et al. (1993) and Bjellqvist et al. (1994); note that the calculation does not implement position-specific adjustments described in the latter paper.
The number of N- and C-terminal groups is taken to be one, unless a value for \code{chains} (number of polypeptide chains) is given in \code{AAcomp}.

Tests for a few proteins (see examples) indicate that \code{GRAVY} and \code{pI} are equal those calculated with the ProtParam tool (\url{https://web.expasy.org/protparam/}; Gasteiger et al., 2005).
}

\references{
Bjellqvist, B., Hughes, G. J., Pasquali, C., Paquet, N., Ravier, F., Sanchez, J.-C., Frutiger, S. and Hochstrasser, D. (1993) The focusing positions of polypeptides in immobilized pH gradients can be predicted from their amino acid sequences. \emph{Electrophoresis} \bold{14}, 1023--1031. \doi{10.1002/elps.11501401163}

Bjellqvist, B. and Basse, B. and Olsen, E. and Celis, J. E. (1994) Reference points for comparisons of two-dimensional maps of proteins from different human cell types defined in a {pH} scale where isoelectric points correlate with polypeptide compositions. \emph{Electrophoresis} \bold{15}, 529--539. \doi{10.1002/elps.1150150171}

Dick, J. M., Yu, M. and Tan, J. (2019) Geochemistry of metagenomes: Redox and salinity gradients and the chemical compositions of proteins. \emph{(preprint)}

Gasteiger, E., Hoogland, C., Gattiker, A., Duvaud, S., Wilkins, M. R., Appel, R. D. and Bairoch, A. (2005) Protein identification and analysis tools on the ExPASy server. In J. M. Walker (Ed.), \emph{The Proteomics Protocols Handbook} (pp. 571--607). Totowa, NJ: Humana Press Inc. \doi{10.1385/1-59259-890-0:571}

Kyte, J. and Doolittle, R. F. (1982) A simple method for displaying the hydropathic character of a protein. \emph{J. Mol. Biol.} \bold{157}, 105--132. \doi{10.1016/0022-2836(82)90515-0}
}

\examples{
# we need CHNOSZ for these examples
require(CHNOSZ)

# for reference, compute ZC of alanine and glycine "by hand"
ZC.Gly <- ZC("C2H5NO2")
ZC.Ala <- ZC("C3H7NO2")
# define the composition of a Gly-Ala-Gly tripeptide
AAcomp <- data.frame(Gly = 2, Ala = 1)
# calculate the ZC of the tripeptide (value: 0.571)
ZC.GAG <- ZCAA(AAcomp)
# this is equal to the carbon-number-weighted average of the amino acids
nC.Gly <- 2 * 2
nC.Ala <- 1 * 3
ZC.average <- (nC.Gly * ZC.Gly + nC.Ala * ZC.Ala) / (nC.Ala + nC.Gly)
stopifnot(all.equal(ZC.GAG, ZC.average))

# compute the per-residue nH2O of Gly-Ala-Gly
basis("QEC")
nH2O.GAG <- species("Gly-Ala-Gly")$H2O
# divide by the length to get residue average (we keep the terminal H-OH)
nH2O.residue <- nH2O.GAG / 3
# compare with the value calculated by H2OAA() (-0.2)
nH2O.H2OAA <- H2OAA(AAcomp, "QEC")
stopifnot(all.equal(nH2O.residue, nH2O.H2OAA))

# calculate GRAVY for lysozyme and ribonuclease
# first get the protein index in CHNOSZ's list of proteins
iprotein <- pinfo(c("LYSC_CHICK", "RNAS1_BOVIN"))
# then get the amino acid compositions
AAcomp <- pinfo(iprotein)
# then calculate GRAVY
Gcalc <- as.numeric(GRAVY(AAcomp))
# these are equal to values obtained with ProtParam on uniprot.org
Gref <- c(-0.472, -0.663)
stopifnot(all.equal(round(Gcalc, 3), Gref))

# calculate pI for a few proteins
iprotein <- pinfo(c("CSG_HALJP", "AMYA_PYRFU", "RNAS1_BOVIN", "LYSC_CHICK"))
AAcomp <- pinfo(iprotein)
pI_calc <- pI(AAcomp)
# reference values calculated with ProtParam on uniprot.org
# CSG_HALJP: residues 35-862 (sequence v1)
# AMYA_PYRFU: residues 2-649 (sequence v2)
# RNAS1_BOVIN: residues 27-150 (sequence v1)
# LYSC_CHICK: residues 19-147 (sequence v1)
pI_ref <- c(3.37, 5.46, 8.64, 9.32)
stopifnot(all.equal(as.numeric(pI_calc), pI_ref))
}
