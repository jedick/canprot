\encoding{UTF-8}
\name{human.aa}
\alias{human.aa}
\title{Get amino acid compositions of human proteins}
\description{
  Get amino acid compositions of human proteins from their UniProt IDs.
}

\usage{
  human.aa(uniprot = NULL, aa_file = NULL,
    stop.if.missing = FALSE, warn.if.duplicated = FALSE)
}

\arguments{
  \item{uniprot}{character, UniProt IDs of proteins}
  \item{aa_file}{character, file name}
  \item{stop.if.missing}{logical, stop with an error if there are UniProt IDs that can't be found?}
  \item{warn.if.duplicated}{logical, emit a warning if duplicate UniProt IDs are detected?}
}

\details{
This function retrieves the amino acid compositions of one or more proteins specified by \code{uniprot}.

This function depends on the amino acid compositions of human proteins, which are stored in the \code{\link{canprot}} environment when the package is attached.
If \code{aa_file} is specified, additional amino acid compositions are read from this file.
This file should be in the same format as \code{\link{human_extra}.csv} in the installation directory of the package.
}

\value{
The function returns a data frame with amino acid compositions of proteins.
}

\examples{
human.aa("P24298")
}

\concept{Amino acid composition}
