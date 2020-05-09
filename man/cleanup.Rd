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
  \item{up2}{logical, TRUE for up-regulated proteins, FALSE for down-regulated proteins}
}

\details{
\code{cleanup} is used in the \code{\link{pdat_}} functions to clean up the dataset given in \code{dat}.
\code{IDcol} is the name of the column that has the UniProt IDs, and \code{up2} indicates the expression change for each protein.
The function removes proteins with unavailable (NA or "") or duplicated IDs.
If \code{up2} is provided, the function also remove unquantified proteins (those that have NA values of \code{up2}) and those with ambiguous expression ratios (up and down for the same ID).
For each operation, a message is printed describing the number of proteins that are \samp{unavailable}, \samp{unquantified}, \samp{ambiguous}, or \samp{duplicated}.

Alternatively, if \code{IDcol} is a logical value, it selects proteins to be unconditionally removed.
}

\seealso{
This function is used extensively by the \code{\link{pdat_}} functions, where it is called after \code{\link{check_IDs}} (if needed).
}

\examples{
# Set up a simple workflow
extdatadir <- system.file("extdata", package="canprot")
datadir <- paste0(extdatadir, "/expression/pancreatic/")
dataset <- "CYD+05"
dat <- read.csv(paste0(datadir, dataset, ".csv.xz"), as.is = TRUE)
up2 <- dat$Ratio..cancer.normal. > 1
# Remove two unavailable and one duplicated proteins
dat <- cleanup(dat, "Entry", up2)
# Now we can calculate the chemical compositions
pcomp <- protcomp(dat$Entry)

# Read another data file
datadir <- paste0(system.file("extdata", package="canprot"), "/expression/colorectal/")
dataset <- "STK+15"
dat <- read.csv(paste0(datadir, "STK+15.csv.xz"), as.is = TRUE)
# Remove unavailable proteins
dat <- cleanup(dat, "uniprot")
# Remove proteins that have less than 2-fold expression ratio
dat <- cleanup(dat, abs(log2(dat$invratio)) < 1)
}

\concept{Protein data}
