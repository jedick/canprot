\encoding{UTF-8}
\name{update_IDs}
\alias{update_IDs}
\title{Update UniProt IDs}
\description{
  Replace old UniProt IDs with current IDs.
}

\usage{
  update_IDs(dat, IDcol, updates_file = NULL)
}

\arguments{
  \item{dat}{data frame, protein expression data}
  \item{IDcol}{character, name of column that has the UniProt IDs}
  \item{updates_file}{character, file name}
}

\details{
The mapping between old and current IDs is taken from the \code{\link{uniprot_updates}} file in the package.
If \code{updates_file} is specified, these ID mappings are included as well.
}

\value{
\code{dat} with updated IDs.
}

\seealso{
\code{\link{cleanup}}
}

\examples{
dat <- data.frame(Entry = "P50224", stringsAsFactors = FALSE)
update_IDs(dat, "Entry")
}

\concept{Protein expression}
