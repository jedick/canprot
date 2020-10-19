\encoding{UTF-8}
\name{check_IDs}
\alias{check_IDs}
\title{Check UniProt IDs}
\description{
Find the first ID for each protein that matches a known UniProt ID.
}

\usage{
  check_IDs(dat, IDcol, aa_file = NULL, updates_file = NULL)
}

\arguments{
  \item{dat}{data frame, protein expression data}
  \item{IDcol}{character, name of column that has the UniProt IDs}
  \item{aa_file}{character, name of file with additional amino acid compositions}
  \item{updates_file}{character, name of file with old to new ID mappings}
}

\details{
\code{check_IDs} is used to check for known UniProt IDs and to update obsolete IDs.
The source IDs should be provided in the \code{IDcol} column of \code{dat}; multiple IDs for one protein can be separated by a semicolon.

The function keeps the first \dQuote{known} ID for each protein, which must be present in one of these groups:
\itemize{
  \item The \code{\link{human_aa}} dataset of amino acid compositions.
  \item Old UniProt IDs that are mapped to new UniProt IDs in \code{\link{uniprot_updates}} or in \code{updates_file} if specified.
  \item IDs of proteins in \code{aa_file}, which lists amino acid compositions in the format described for \code{\link{human_aa}} (see \code{extdata/protein/human_extra.csv} for an example and \code{\link[CHNOSZ]{thermo}$protein} for more details).
}
}

\value{
\code{dat} is returned with possibly changed values in the column designated by \code{IDcol}; old IDs are replaced with new ones, the first known ID for each protein is kept, then proteins with no known IDs are assigned \code{NA}.
}

\seealso{
This function is used by the \code{\link{pdat_}} functions, where it is called before \code{\link{cleanup}}.
}

\examples{
# Make up some data for this example
ID <- c("P61247;PXXXXX", "PYYYYY;P46777;P60174", "PZZZZZ")
dat <- data.frame(ID = ID, stringsAsFactors = FALSE)
# Get the first known ID for each protein; the third one is NA
check_IDs(dat, "ID")

# Update an old ID
dat <- data.frame(Entry = "P50224", stringsAsFactors = FALSE)
check_IDs(dat, "Entry")
}

\concept{Protein data}
