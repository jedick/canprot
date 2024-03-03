\encoding{UTF-8}
\name{metrics}
\alias{metrics}
\alias{nH2O}
\alias{nO2}
\alias{Zc}
\alias{GRAVY}
\alias{pI}
\alias{MW}
\alias{V0}
\alias{pV0}
\alias{plength}
\alias{HC}
\alias{NC}
\alias{OC}
\alias{SC}
\alias{cplab}
\title{Calculate Chemical Metrics for Proteins}
\description{
Calculate chemical metrics for proteins from their amino acid compositions.
}

\usage{
  Zc(AAcomp, ...)
  nH2O(AAcomp, basis = "QEC", terminal_H2O = 0)
  nO2(AAcomp, basis = "QEC", ...)
  GRAVY(AAcomp, ...)
  pI(AAcomp, terminal_H2O = 1, ...)
  MW(AAcomp, terminal_H2O = 0, ...)
  V0(AAcomp, terminal_H2O = 0, ...)
  pV0(AAcomp, terminal_H2O = 0, ...)
  HC(AAcomp, ...)
  NC(AAcomp, ...)
  OC(AAcomp, ...)
  SC(AAcomp, ...)
  plength(AAcomp, ...)
  cplab
}

\arguments{
  \item{AAcomp}{data frame, amino acid compositions}
  \item{...}{ignored additional arguments}
  \item{basis}{character, set of basis species}
  \item{terminal_H2O}{numeric, number of pairs of terminal groups}
}

\details{
Columns in \code{AAcomp} should be named with the three-letter abbreviations for the amino acids (\samp{Ala}, \samp{Arg}, \ldots).
Case-insensitive matching matching of the abbreviations is used (e.g. \samp{ALA} is the same as \samp{ala}).

Metrics are normalized per amino acid residue except for \code{Zc}, \code{GRAVY}, \code{pI}, and \code{plength}.
The contribution of protein terminal groups (-H and -OH) to relevant metrics is turned off by default, except for isoelectric point \code{pI}.
Set \code{terminal_H2O} to 1 (or to the number of polypeptide chains, if greater than one) to include their contribution.

The metrics are described below:

\describe{
\item{\code{Zc}}{Average oxidation state of carbon (\Zc) (Dick, 2014).
  This metric is independent of the choice of basis species.
  Note that \Zc is normalized by number of carbon atoms, not by number of residues.
}
\item{\code{nH2O}}{Stoichiometric hydration state (\nH2O) per residue.
  The available \code{basis} species are:
  \itemize{
    \item{\samp{QEC} - glutamine, glutamic acid, cysteine, \H2O, \O2 (Dick et al., 2020)}
    \item{\samp{QCa} - glutamine, cysteine, acetic acid, \H2O, \O2}
  }
}
\item{\code{nO2}}{Stoichiometric oxidation state (\nO2) per residue.
  The basis species also affect this calculation.
}
\item{\code{GRAVY}}{Grand average of hydropathy.
  Values of the hydropathy index for individual amino acids are from Kyte and Doolittle (1982).
}
\item{\code{pI}}{Isoelectric point.
  The net charge for each ionizable group was pre-calculated from pH 0 to 14 at intervals of 0.01.
  The isoelectric point is found as the pH where the sum of charges of all groups in the protein is closest to zero.
  The \pK values for the terminal groups and sidechains are taken from Bjellqvist et al. (1993) and Bjellqvist et al. (1994); note that the calculation does not implement position-specific adjustments described in the latter paper.
  The number of N- and C-terminal groups is taken from \code{terminal_H2O}.
}
\item{\code{MW}}{Molecular weight per residue.}
\item{\code{V0}}{Standard molal volume.
  The values are derived from group contributions of amino acid sidechains and protein backbones (Dick et al., 2006).
}
\item{\code{plength}}{Protein length (number of amino acid residues).}
\item{\code{HC}}{H/C ratio (not counting terminal -H and -OH groups).}
\item{\code{NC}}{N/C ratio.}
\item{\code{OC}}{O/C ratio (not counting terminal -H and -OH groups).}
\item{\code{SC}}{S/C ratio.}
}

\code{...} is provided to permit \code{\link{get}} or \code{\link{do.call}} constructions with the same arguments for all metrics.
For instance, a \code{terminal_H2O} argument can be suppled to either \code{Zc} or \code{nH2O}, but it only has an effect on the latter.

\code{cplab} is a list of formatted labels for each of the chemical metrics listed here.
}

\references{
Bjellqvist B, Hughes GJ, Pasquali C, Paquet N, Ravier F, Sanchez J-C, Frutiger S, Hochstrasser D. 1993. The focusing positions of polypeptides in immobilized pH gradients can be predicted from their amino acid sequences. \emph{Electrophoresis} \bold{14}: 1023--1031. \doi{10.1002/elps.11501401163}

Bjellqvist B, Basse B, Olsen E, Celis JE. 1994. Reference points for comparisons of two-dimensional maps of proteins from different human cell types defined in a {pH} scale where isoelectric points correlate with polypeptide compositions. \emph{Electrophoresis} \bold{15}: 529--539. \doi{10.1002/elps.1150150171}

Dick JM, LaRowe DE, Helgeson HC. 2006. Temperature, pressure, and electrochemical constraints on protein speciation: Group additivity calculation of the standard molal thermodynamic properties of ionized unfolded proteins. \emph{Biogeosciences} \bold{3}(3): 311--336. \doi{10.5194/bg-3-311-2006}

Dick JM. 2014. Average oxidation state of carbon in proteins. \emph{J. R. Soc. Interface} \bold{11}: 20131095. \doi{10.1098/rsif.2013.1095}

Dick JM, Yu M, Tan J. 2020. Uncovering chemical signatures of salinity gradients through compositional analysis of protein sequences. \emph{Biogeosciences} \bold{17}: 6145--6162. \doi{10.5194/bg-17-6145-2020}

Kyte J, Doolittle RF. 1982. A simple method for displaying the hydropathic character of a protein. \emph{J. Mol. Biol.} \bold{157}: 105--132. \doi{10.1016/0022-2836(82)90515-0}

}

\seealso{
For calculation of \Zc from an elemental formula (instead of amino acid composition), see the \code{\link[CHNOSZ]{ZC}} function in \pkg{CHNOSZ}.
\code{\link{calc.metrics}} is a wrapper to calculate one or more metrics specified in an argument.
}

\examples{
# Amino acid composition of a tripeptide (Gly-Ala-Gly)
aa <- data.frame(Ala = 1, Gly = 2)
# Calculate Zc, nH2O, and length
Zc(aa)
nH2O(aa)
plength(aa)

# Make a plot with formatted labels
plot(Zc(aa), nH2O(aa), xlab = cplab$Zc, ylab = cplab$nH2O)
}
