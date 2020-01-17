\encoding{UTF-8}
\name{cleanup}
\alias{cleanup}
\title{Clean up Data}
\description{
  Remove proteins with unavailable IDs, ambiguous expression ratios, and duplicated IDs.
}

\usage{
  cleanup(dat, IDcol, up2 = NULL)
}

\arguments{
  \item{dat}{data frame, protein expression data}
  \item{IDcol}{character, name of column that has the UniProt IDs}
  \item{up2}{logical, TRUE for proteins that are up-expressed}
}

\details{
\code{cleanup} is used in the \code{\link{pdat_}} functions to \dQuote{clean up} a dataset in \code{dat}.
\code{IDcol} should identify the column that has the UniProt IDs, and \code{up2} indicates the expression change for proteins (TRUE for up-expressed).
The function removes proteins with unavailable or duplicated IDs.
If \code{up2} is provided, the function also remove unquantified proteins (those that have NA values of \code{up2}) and those with ambiguous expression ratios (up and down for the same ID).
For each operation, a message is printed describing the number of proteins that are \samp{unavailable}, \samp{unquantified}, \samp{ambiguous}, and \samp{duplicated}.

Alternatively, if \code{IDcol} is a logical value, it selects proteins to be unconditionally removed.
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
dat <- cleanup(dat, "Entry", up2)
pcomp <- protcomp(dat$Entry)

# using cleanup() to remove "manually" selected proteins
datadir <- paste0(system.file("extdata", package="canprot"), "/expression/colorectal/")
dataset <- "STK+15"
dat <- read.csv(paste0(datadir, "STK+15.csv.xz"), as.is=TRUE)
dat <- cleanup(dat, is.na(dat$uniprot))
}

\concept{Protein expression}
