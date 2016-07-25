\encoding{UTF-8}
\name{remove_entries}
\alias{remove_entries}
\title{Remove Data Entries}
\description{
  Remove rows of a data table matching a condition, and print an informative message.
}

\usage{
  remove_entries(dat, irm, dataset, description) 
}

\arguments{
  \item{dat}{data frame}
  \item{irm}{logical, which rows to remove}
  \item{dataset}{character, name of dataset used in message}
  \item{description}{character, description used in message}
}

\details{
This function is used by the \code{\link{pdat_}} functions to \dQuote{clean up} particular datasets.
Examples of the \code{description} in \code{\link{pdat_CRC}} are \samp{missing}, \samp{conflicting}, and \samp{duplicated}.
This function is made accessible to the user who may wish to write their own \code{\link{pdat_}} functions.
}

\examples{
extdatadir <- paste0(system.file("extdata", package="canprot"), "/")
dataset <- "STK+15"
# an excerpt from pdat_CRC() for dataset="STK+15"
dat <- read.csv(paste0(extdatadir, "STK+15.csv"), as.is=TRUE)
dat <- remove_entries(dat, is.na(dat$uniprot), dataset, "missing")
# the above prints "STK+15: dropping 3 missing proteins"
}

\concept{Protein expression}
