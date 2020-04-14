\encoding{UTF-8}
\name{recomp}
\alias{recomp}
\title{Recompute Protein Compositions}
\description{
  Recompute compositional metrics given the existing output of one of the \code{\link{pdat_}} functions.
}

\usage{
  recomp(pdat, basis = "rQEC")
}

\arguments{
  \item{pdat}{list, result from a \code{\link{pdat_}} function}
  \item{basis}{character, keyword for basis species to use}
}

\details{
This function recomputes a pdat object for a different \samp{basis}.
This avoids the need to regenerate the list of UniProt IDs and search for amino acid compositions, and is used to speed up the processing of the vignettes.
}

\examples{
# process the first dataset for prostate cancer
datasets <- pdat_prostate()
pdat <- pdat_prostate(datasets[1], basis = "QEC")
# recompute the compositions using amino acid biosynthetic reactions
pdat2 <- recomp(pdat, "biosynth")
# compare the results
nO2_QEC <- pdat$pcomp$residue.basis[, "O2"]
nO2_biosynth <- pdat2$pcomp$residue.basis[, "O2"]
plot(nO2_QEC, nO2_biosynth)
# the slower way (not that it makes much difference for this example)
pdat2.slow <- pdat_prostate(datasets[1], basis = "biosynth")
stopifnot(identical(pdat2, pdat2.slow))
}

\concept{Chemical composition}
