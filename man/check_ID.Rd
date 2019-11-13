\encoding{UTF-8}
\name{check_ID}
\alias{check_ID}
\title{Check UniProt IDs}
\description{
Do the IDs correspond to known UniProt IDs?
}

\usage{
  check_ID(ID, aa_file = NULL, updates_file = NULL)
}

\arguments{
  \item{ID}{character or list, candidate UniProt IDs}
  \item{aa_file}{character, file name}
  \item{updates_file}{character, file name}
}

\details{
\code{check_ID} is used to check candidate \code{ID}s supplied in a character or list object.
Multiple IDs can be separated by a semicolon.
\dQuote{Known} IDs are those that are present in the \code{\link{human_aa}} dataset of amino acid compositions.
Known IDs also include old UniProt IDs that are mapped to new UniProt IDs in \code{\link{uniprot_updates}} and in \code{updates_file} if specified.

If \code{aa_file} is specified, additional amino acid compositions to be considered are read from this file, which should be in the same format as e.g. \code{\link{human_extra}.csv} (see also \code{\link[CHNOSZ]{thermo}$protein}).
}

\value{
The function returns the IDs in a list (dropping the semicolons, but reflecting the original arrangement), replacing any unknown IDs with \code{NA}.
}

\examples{
# the function replaces the 3 non-existent IDs with NA
check_ID(c("P61247;PXXXXX", "PYYYYY;P46777;P60174", "PZZZZZ"))
}

\concept{Amino acid composition}
