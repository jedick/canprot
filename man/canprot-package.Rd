\name{canprot-package}
\alias{canprot-package}
\docType{package}
\title{Differential Expression of Proteins in Cancer}
\description{
  \pkg{canprot} is a package for exploration of compositional data of proteomes and thermodynamic analysis of proteomic transformations using the concepts of chemical components (basis species) and chemical affinity (negative of Gibbs energy).
}

\section{Overview}{
  This package includes datasets for differential expression of proteins in colorectal cancer, pancreatic cancer, hypoxia, and hyperosmotic stress.
  All datasets use UniProt IDs, which have been added if not present in the original publications. 
  The sources of data are listed in the vignettes and are further described in Dick (2016) and (2017).
  The functions in this package were derived from the code included as supplemental information for those studies.
}

\examples{
# list all of the data files for protein expression
exprdata <- system.file("extdata/expression", package="canprot")
exprfiles <- dir(exprdata, recursive=TRUE)
print(exprfiles)
# get the reference keys from the filenames
refkeys <- gsub(".csv", "", sapply(strsplit(exprfiles, "/"), "[", 2))
# find the reference keys in the UniProt updates file
uu <- get("uniprot_updates", canprot)
update_keys <- unique(unlist(strsplit(uu$source, ";")))
# find the reference keys in the extra human amino acid composition file
he <- get("human_extra", canprot)
extra_keys <- unique(unlist(strsplit(he$ref, ";")))
# list the unused keys (these should be empty when the package is released)
setdiff(update_keys, refkeys)
setdiff(extra_keys, refkeys)
}

\references{
  Dick, Jeffrey M. (2016) Proteomic indicators of oxidation and hydration state in colorectal cancer. \emph{PeerJ} \bold{4}, e2238. \doi{10.7717/peerj.2238}

  Dick, Jeffrey M. (2017) Chemical composition and the potential for proteomic transformation in cancer, hypoxia, and hyperosmotic stress. \emph{PeerJ} \bold{5}, e3421. \doi{10.7717/peerj.3421}
}

\keyword{package}
