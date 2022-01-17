\encoding{UTF-8}
\name{protcomp}
\alias{protcomp}
\title{Amino Acid Compositions}
\description{
  Get amino acid compositions of proteins.
}

\usage{
  protcomp(uniprot = NULL, aa = NULL, aa_file = NULL)
}

\arguments{
  \item{uniprot}{character, UniProt IDs of proteins}
  \item{aa}{data frame, amino acid compositions}
  \item{aa_file}{character, file name}
}

\details{
This function retrieves the amino acid compositions of one or more proteins specified by \code{uniprot}.

This function depends on the amino acid compositions of human proteins, which are stored in the \code{\link{human}} environment when the package is attached.
If \code{aa_file} is specified, additional amino acid compositions to be considered are read from this file, which should be in the same format as e.g. \code{\link{human_extra}.csv} (see also \code{\link[CHNOSZ]{thermo}$protein}).
Alternatively, the amino acid compositions can be given in \code{aa}, bypassing the search step.
}

\value{
The function returns a list with elements \code{uniprot} (UniProt IDs as given in the arguments) and \code{aa} (amino acid compositions of the proteins).
}

\seealso{
\code{\link{cleanup}}
}

\examples{
protcomp("P24298")
}

\concept{Amino acid composition}
