\encoding{UTF-8}
\name{cleanup}
\alias{cleanup}
\title{Clean up Data}
\description{
  Remove proteins with unavailable IDs, ambiguous expression ratios, and duplicated IDs.
}

\usage{
  cleanup(dat, IDcol, dataset, up2)
}

\arguments{
  \item{dat}{data frame, protein expression data}
  \item{IDcol}{character, name of column that has the UniProt IDs}
  \item{dataset}{character, name of dataset}
  \item{up2}{logical, which proteins are up-regulated}
}

\value{
\code{dat} excluding any proteins with unavailable or duplicated IDs, or ambiguous expression ratios (up and down for the same ID).
}

\seealso{
\code{\link{update_IDs}}
}

\examples{
# cleanup() used in a data analysis workflow
# here, it removes 2 unavailable and 1 duplicated protein
extdatadir <- system.file("extdata", package="canprot")
datadir <- paste0(extdatadir, "/expression/pancreatic/")
dataset <- "CYD+05"
dat <- read.csv(paste0(datadir, dataset, ".csv.xz"), as.is=TRUE)
up2 <- dat$Ratio..cancer.normal. > 1
dat <- cleanup(dat, "Entry", dataset, up2)
pcomp <- protcomp(dat$Entry)
}

\concept{Protein expression}
