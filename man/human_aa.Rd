\encoding{UTF-8}
\name{human_aa}
\alias{human_aa}
\title{Get amino acid compositions of human proteins}
\description{
  Get amino acid compositions of human proteins from their UniProt IDs.
}

\usage{
  human_aa(uniprot = NULL, aa_file = NULL,
    stop_if_missing = FALSE, warn_if_duplicated = FALSE)
}

\arguments{
  \item{uniprot}{character, UniProt IDs of proteins}
  \item{aa_file}{character, file name}
  \item{stop_if_missing}{logical, stop with an error if there are UniProt IDs that can't be found?}
  \item{warn_if_duplicated}{logical, emit a warning if duplicate UniProt IDs are detected?}
}

\details{
This function retrieves the amino acid compositions of one or more proteins specified by \code{uniprot}.

This function depends on the amino acid compositions of human proteins, which are stored in the \code{\link{canprot}} environment when the package is attached.
If \code{aa_file} is specified, additional amino acid compositions are read from this file.
This file should be in the same format as \code{\link{human.extra}.csv} in the installation directory of the package.
}

\value{
The function returns a data frame with amino acid compositions of proteins.
}

\examples{
human_aa("P24298")
}

\concept{Amino acid composition}
