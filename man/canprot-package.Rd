\name{canprot-package}
\alias{canprot-package}
\docType{package}
\title{Differential Expression of Proteins in Cancer}
\description{
  \pkg{canprot} is a package for exploration of compositional data of proteomes and thermodynamic analysis of proteomic transformations using the concepts of chemical components (basis species) and chemical affinity (negative of Gibbs energy).
}

\section{Overview}{
  Datasets of differential expression of proteins are included for colorectal cancer (CRC), pancreatic cancer, and hypoxia.
  Please see the vignettes, starting with \href{../doc/CRC.html}{CRC.Rmd}, for the data sources and examples of compositional and potential diagrams.
  This package was initiated using code and datasets for CRC that were assembled for the study of Dick (2016).
}

\section{Acknowledgements}{
  I thank Dr. Ming-Chih Lai for provided unpublished data related to the study of Lai et al. (2016) and giving permission to distribute it.
}

\examples{
# list all of the data files for protein expression
exprdata <- system.file("extdata/expression", package="canprot")
exprfiles <- dir(exprdata, recursive=TRUE)
print(exprfiles)
# get the reference keys from the filenames
refkeys <- gsub(".csv", "", sapply(strsplit(exprfiles, "/"), "[", 2))
# find the reference keys in the UniProt updates file
data(canprot)
update_keys <- unique(unlist(strsplit(uniprot_updates$source, ";")))
# find the reference keys in the extra human amino acid composition file
extra_keys <- unique(unlist(strsplit(human_extra$ref, ";")))
# list the unused keys (these should be cleaned up before the package is released)
setdiff(update_keys, refkeys)
setdiff(extra_keys, refkeys)
}

\references{
  Dick, Jeffrey M. (2016) Proteomic indicators of oxidation and hydration state in colorectal cancer. \emph{PeerJ} \bold{4}, e2238. \doi{10.7717/peerj.2238}

  Lai, Ming-Chih and Chang, Chiao-May and Sun, H. Sunny (2016)  Hypoxia induces autophagy through translational up-regulation of lysosomal proteins in human colon cancer cells. \emph{PLoS ONE} \bold{11}(4), e0153627. \doi{10.1371/journal.pone.0153627}
}


\keyword{package}
