\encoding{UTF-8}
\name{PS}
\alias{PS}
\title{Retrieve phylostrata for given UniProt IDs}
\description{
This retrieves the phylostrata according to Trigos et al. (2017)
}

\usage{
  PS(uniprot)
}

\arguments{
  \item{uniprot}{character, UniProt accession numbers}
}

\details{
The data file at \code{extdata/phylostrata/TPPG17.csv} has these columns: \samp{GeneID} (gene name), \samp{Entrez}, \samp{Entry}, \samp{Phylostrata}.
Except for \samp{Entry}, the values are taken from Dataset S1 of Trigos et al. (2017).
UniProt acession numbers in \samp{Entry} were generated using the UniProt mapping tool first for \samp{Entrez}, followed by \samp{GeneID} for the unmatched genes.
\samp{Entry} is NA for genes that remain unmatched to any proteins after both mapping steps.
}

\references{
Trigos, A. S. and Pearson, R. B. and Papenfuss, A. T. and Goode, D. L. (2017) Altered interactions between unicellular and multicellular genes drive hallmarks of transformation in a diverse range of solid tumors. \emph{Proc. Natl. Acad. Sci.} \bold{114}, 6406--6411. \doi{10.1073/pnas.1617743114}
}

\seealso{
\code{\link{metrics}}
}

\examples{
# plot nH2O vs phylostrata for some proteins
pd <- get_pdat("JKMF10")
get_comptab(pd, "PS", plot.it = TRUE)
}
