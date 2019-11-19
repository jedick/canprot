\encoding{UTF-8}
\name{check_IDs}
\alias{check_IDs}
\title{Check UniProt IDs}
\description{
Do the IDs correspond to known UniProt IDs?
}

\usage{
  check_IDs(dat, IDcol, aa_file = NULL, updates_file = NULL)
}

\arguments{
  \item{dat}{data frame, protein expression data}
  \item{IDcol}{character, name of column that has the UniProt IDs}
  \item{aa_file}{character, file name}
  \item{updates_file}{character, file name}
}

\details{
\code{check_IDs} is used to check candidate UniProt IDs in the \code{IDcol} column of \code{dat}.
Multiple IDs for one protein can be separated by a semicolon.
\dQuote{Known} IDs are those that are present in the \code{\link{human_aa}} dataset of amino acid compositions.
Known IDs also include old UniProt IDs that are mapped to new UniProt IDs in \code{\link{uniprot_updates}} and in \code{updates_file} if specified.

If \code{aa_file} is specified, additional amino acid compositions to be considered are read from this file, which should be in the same format as e.g. \code{\link{human_extra}.csv} (see also \code{\link[CHNOSZ]{thermo}$protein}).
}

\value{
The function returns \code{dat} with possibly updated values in \code{IDcol}.
All unknown IDs are replaced with \code{NA}.
}

\examples{
ID <- c("P61247;PXXXXX", "PYYYYY;P46777;P60174", "PZZZZZ")
dat <- data.frame(ID = ID, stringsAsFactors = FALSE)
check_IDs(dat, "ID")
}

\concept{Amino acid composition}
