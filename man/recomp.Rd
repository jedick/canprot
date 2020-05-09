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
This avoids the need to regenerate the list of UniProt IDs and search for amino acid compositions.
It was previously used to speed up the processing of the vignettes.
}

\examples{
# process the first dataset for prostate cancer using the QEC basis species
datasets <- pdat_prostate()
pdat <- pdat_prostate(datasets[1], basis = "QEC")
# recompute the compositions using the rQEC values
pdat2 <- recomp(pdat, "rQEC")
# compare the results
nH2O_QEC <- pdat$pcomp$residue.basis[, "H2O"]
nH2O_rQEC <- pdat2$pcomp$residue.basis[, "H2O"]
plot(nH2O_QEC, nH2O_rQEC)
# the slower way (not that it makes much difference for this example)
pdat2.slow <- pdat_prostate(datasets[1], basis = "rQEC")
stopifnot(identical(pdat2, pdat2.slow))
}

\concept{Chemical composition}
