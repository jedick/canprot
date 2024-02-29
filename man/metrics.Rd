\encoding{UTF-8}
\name{metrics}
\alias{metrics}
\alias{nH2O}
\alias{nO2}
\alias{Zc}
\alias{GRAVY}
\alias{pI}
\alias{MW}
\alias{plength}
\alias{basis.text}
\alias{cplab}
\title{Calculate Chemical Metrics for Proteins}
\description{
These functions calculate chemical metrics of proteins given a data frame of amino acid compositions.
}

\usage{
  Zc(AAcomp, ...)
  nH2O(AAcomp, basis = getOption("basis"), terminal_H2O = 0)
  nO2(AAcomp, basis = getOption("basis"), ...)
  GRAVY(AAcomp, ...)
  pI(AAcomp, ...)
  MW(AAcomp, ...)
  plength(AAcomp, ...)
  basis.text(basis)
  cplab
}

\arguments{
  \item{AAcomp}{data frame, amino acid compositions}
  \item{...}{ignored additional arguments}
  \item{basis}{character, basis species}
  \item{terminal_H2O}{numeric, number of pairs of terminal groups}
}

\details{
Columns in \code{AAcomp} should be named with the three-letter abbreviations for the amino acids (\samp{Ala}, \samp{Arg}, \ldots).
Abbreviations are matched without regard to case (e.g. \samp{ALA} is the same as \samp{ala}).

\code{...} is provided to permit \code{\link{get}} or \code{\link{do.call}} constructions with the same arguments for all metrics.
For instance, a \code{terminal_H2O} argument can be suppled to either \code{Zc} or \code{nH2O}, but it only has an effect on the latter.

The metrics are described below:

\describe{
\item{\code{Zc}}{Average oxidation state of carbon (\Zc) (Dick, 2014).
  This metric is independent of the choice of basis species.
  Note that \Zc is normalized by number of carbon atoms, but \nH2O and \nO2 are normalized by number of residues.
}
\item{\code{nH2O}}{Stoichiometric hydration state (\nH2O) per residue.
  The available \code{basis} species are:
  \itemize{
    \item{\samp{QEC} - glutamine, glutamic acid, cysteine, \H2O, \O2 (Dick et al., 2020) (this is the default for \code{getOption("basis")})}
    \item{\samp{QCa} - glutamine, cysteine, acetic acid, \H2O, \O2}
    \item{Any other valid basis specification for \code{\link[CHNOSZ]{basis}}, such as \samp{CHNOS} for \CO2, \NH3, \H2S, \H2O, and \O2}
  }
  The default of 0 for \code{terminal_H2O} means to ignore the terminal -H and -OH groups of proteins.
  This can be changed to 1 (or the number of polypeptide chains) to add one (or more) \nH2O before normalizing by the number of residues.
}
\item{\code{nO2}}{Stoichiometric oxidation state (\nO2) per residue.
  The basis species also affect this calculation.
}
\item{\code{GRAVY}}{Grand average of hydropathicity.
  Values of the hydropathy index for individual amino acids are from Kyte and Doolittle (1982).
}
\item{\code{pI}}{Isoelectric point.
  The net charge for each ionizable group was pre-calculated from pH 0 to 14 at intervals of 0.01.
  The isoelectric point is found as the pH where the sum of charges of all groups in the protein is closest to zero.
  The \pK values for the terminal groups and sidechains are taken from Bjellqvist et al. (1993) and Bjellqvist et al. (1994); note that the calculation does not implement position-specific adjustments described in the latter paper.
  The number of N- and C-terminal groups is taken to be one, unless a value for \code{chains} (number of polypeptide chains) is given in \code{AAcomp}.
}
\item{\code{MW}}{Molecular weight per residue (not including terminal groups).}
\item{\code{plength}}{Protein length (number of amino acid residues).}
}

Tests for a few proteins (see examples) indicate that \code{GRAVY} and \code{pI} are equal those calculated with the ProtParam tool (\url{https://web.expasy.org/protparam/}; Gasteiger et al., 2005).

\code{basis.text} is used in the vignettes to generate a textual description of the names of the basis species, except \H2O and \O2, for one of the keywords \samp{QEC} or \samp{QCa}.

\code{cplab} is a list of formatted labels for chemical metrics listed in \code{\link{metrics}}.
It is an exported object, available to the user and other packages.
}

\references{
Bjellqvist, B., Hughes, G. J., Pasquali, C., Paquet, N., Ravier, F., Sanchez, J.-C., Frutiger, S. and Hochstrasser, D. (1993) The focusing positions of polypeptides in immobilized pH gradients can be predicted from their amino acid sequences. \emph{Electrophoresis} \bold{14}, 1023--1031. \doi{10.1002/elps.11501401163}

Bjellqvist, B. and Basse, B. and Olsen, E. and Celis, J. E. (1994) Reference points for comparisons of two-dimensional maps of proteins from different human cell types defined in a {pH} scale where isoelectric points correlate with polypeptide compositions. \emph{Electrophoresis} \bold{15}, 529--539. \doi{10.1002/elps.1150150171}

Dick, J. M. (2014) Average oxidation state of carbon in proteins. \emph{J. R. Soc. Interface} \bold{11}, 20131095. \doi{10.1098/rsif.2013.1095}

Dick, J. M., Yu, M. and Tan, J. (2020) Uncovering chemical signatures of salinity gradients through compositional analysis of protein sequences. \emph{Biogeosciences} \bold{17}, 6145--6162. \doi{10.5194/bg-17-6145-2020}

Gasteiger, E., Hoogland, C., Gattiker, A., Duvaud, S., Wilkins, M. R., Appel, R. D. and Bairoch, A. (2005) Protein identification and analysis tools on the ExPASy server. In J. M. Walker (Ed.), \emph{The Proteomics Protocols Handbook} (pp. 571--607). Totowa, NJ: Humana Press Inc. \doi{10.1385/1-59259-890-0:571}

Kyte, J. and Doolittle, R. F. (1982) A simple method for displaying the hydropathic character of a protein. \emph{J. Mol. Biol.} \bold{157}, 105--132. \doi{10.1016/0022-2836(82)90515-0}
}

\seealso{
For calculation of \Zc from an elemental formula (instead of amino acid composition), see the \code{\link[CHNOSZ]{ZC}} function in \pkg{CHNOSZ}.
}

\examples{
# Amino acid composition of a tripeptide (Gly-Ala-Gly)
aa <- data.frame(Ala = 1, Gly = 2)
# Calculate Zc, nH2O, and length
Zc(aa)
nH2O(aa)
plength(aa)
}

\concept{Chemical metrics}
