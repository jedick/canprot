\encoding{UTF-8}
\name{PS}
\alias{PS}
\title{Retrieve phylostrata for given UniProt IDs}
\description{
This retrieves the phylostrata for protein-coding genes according to Liebeskind et al. (2016) or Trigos et al. (2017).
}

\usage{
  PS(uniprot, source = "TPPG17")
}

\arguments{
  \item{uniprot}{character, UniProt accession numbers}
  \item{source}{character, \samp{TPPG17} or \samp{LMM16}}
}

\details{
The data file at \code{extdata/phylostrata/TPPG17.csv.xz} has these columns: \samp{GeneID} (gene name), \samp{Entrez}, \samp{Entry}, \samp{Phylostrata}.
Except for \samp{Entry}, the values are from Dataset S1 of Trigos et al. (2017).
UniProt acession numbers in \samp{Entry} were generated using the UniProt mapping tool first for \samp{Entrez}, followed by \samp{GeneID} for the unmatched genes.
\samp{Entry} is NA for genes that remain unmatched to any proteins after both mapping steps.

The data file at \code{extdata/phylostrata/LMM16.csv.xz} has these columns: \samp{UniProt}, \samp{modeAge}, \samp{PS}.
The data are from file \code{main_HUMAN.csv} in Gene-Ages v1.0 (\url{https://zenodo.org/record/51708}; Liebeskind et al. (2016)).
The modeAges were converted to phylostrata values 1-8 (\samp{PS} column) in this order: Cellular_organisms, Euk_Archaea, Euk+Bac, Eukaryota, Opisthokonta, Eumetazoa, Vertebrata, Mammalia.

Specifying \samp{PS} as a variable name in \code{\link{get_comptab}} retrieve the phylostrata using this function.
}

\references{
Trigos, A. S. and Pearson, R. B. and Papenfuss, A. T. and Goode, D. L. (2017) Altered interactions between unicellular and multicellular genes drive hallmarks of transformation in a diverse range of solid tumors. \emph{Proc. Natl. Acad. Sci.} \bold{114}, 6406--6411. \doi{10.1073/pnas.1617743114}

Liebeskind, B. J. and McWhite, Claire J. and Marcotte, E. M. (2016) Towards consensus gene ages. \emph{Genome Biol. Evol.} \bold{8}, 1812--1823. \doi{10.1093/gbe/evw113}
}

\seealso{
\code{\link{metrics}}
}

\examples{
# plot nH2O vs phylostrata
pd <- pdat_CRC("JKMF10")
# (Trigos et al. data)
get_comptab(pd, "PS", plot.it = TRUE)
# (Liebeskind et al. data)
get_comptab(pd, "PS", plot.it = TRUE, PS_source = "LMM16")

# compare the two sources
TPPG17 <- read.csv(system.file("extdata/phylostrata/TPPG17.csv.xz", package = "canprot"))
LMM16 <- read.csv(system.file("extdata/phylostrata/LMM16.csv.xz", package = "canprot"))
IDs <- intersect(TPPG17$Entry, LMM16$UniProt)
PS_TPPG17 <- TPPG17$Phylostrata[match(IDs, TPPG17$Entry)]
PS_LMM16 <- LMM16$PS[match(IDs, LMM16$UniProt)]
plot(jitter(PS_TPPG17), jitter(PS_LMM16), pch = ".")
}
