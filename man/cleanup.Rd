\encoding{UTF-8}
\name{cleanup}
\alias{cleanup}
\alias{remove_entries}
\title{Clean up Data}
\description{
  Remove proteins with unavailable IDs, ambiguous expression ratios, and duplicated IDs.
}

\usage{
  cleanup(dat, IDcol, dataset, up2)
  remove_entries(dat, irm, dataset, description) 
}

\arguments{
  \item{dat}{data frame, protein expression data}
  \item{IDcol}{character, name of column that has the UniProt IDs}
  \item{dataset}{character, name of dataset}
  \item{up2}{logical, which proteins are up-regulated}
  \item{irm}{logical, which rows to remove}
  \item{description}{character, description used in message}
}

\details{
\code{cleanup} is used in the \code{\link{pdat_}} functions to \dQuote{clean up} particular datasets.

\code{remove_entries} is the \dQuote{low-level} function that actually removes the rows, and prints a message.
The values of \code{description} that are used in \code{cleanup} are \samp{unavailable}, \samp{unquantified}, \samp{ambiguous}, and \samp{duplicated}.
}

\value{
\code{dat} excluding any unquantified proteins (those that have NA values of up2) and any proteins with unavailable or duplicated IDs or ambiguous expression ratios (up and down for the same ID).
}

\seealso{
\code{\link{check_IDs}} (should be used before \code{cleanup}).
}

\examples{
# cleanup() used in a simple workflow
# here, it removes 2 unavailable and 1 duplicated protein
extdatadir <- system.file("extdata", package="canprot")
datadir <- paste0(extdatadir, "/expression/pancreatic/")
dataset <- "CYD+05"
dat <- read.csv(paste0(datadir, dataset, ".csv.xz"), as.is=TRUE)
up2 <- dat$Ratio..cancer.normal. > 1
dat <- cleanup(dat, "Entry", dataset, up2)
pcomp <- protcomp(dat$Entry)

# using remove_entries() "by hand"
datadir <- paste0(system.file("extdata", package="canprot"), "/expression/CRC/")
dataset <- "STK+15"
dat <- read.csv(paste0(datadir, "STK+15.csv.xz"), as.is=TRUE)
dat <- remove_entries(dat, is.na(dat$uniprot), dataset, "unavailable")
}

\concept{Protein expression}
