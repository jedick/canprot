\name{canprot-package}
\alias{canprot-package}
\docType{package}
\title{Differential Expression of Proteins in Cancer}
\description{
  \pkg{canprot} is a package for analysis of the chemical compositions of proteins from their amino acid compositions.
  The package compiles datasets for differentially expressed proteins in cancer and cell culture conditions from over 200 studies.
  It also integrates thermodynamic analysis of differentially expressed proteins using the \pkg{CHNOSZ} package.
}

\section{Overview}{
  This package includes datasets for differential expression of proteins in six cancer types (breast, colorectal, liver, lung, pancreatic, prostate), and four cell culture conditions (hypoxia, hyperosmotic stress, secreted proteins in hypoxia, and 3D compared to 2D growth conditions).
  The vignettes have plots for each cancer type and cell culture condition and references for all data sources used.
  All datasets use UniProt IDs; if not present in the original publications they have been added using the UniProt mapping tool (\url{https://www.uniprot.org/mapping/}).

  The functions in this package were originally derived from the code provided in the supplemental information for Dick (2016) and (2017).
  Updated data compilations and revised plotting functions and vignettes were developed for the studies of Dick et al. (2020) and Dick (2020).
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

\seealso{
  The \pkg{JMDplots} package on GitHub (\url{https://github.com/jedick/JMDplots}) has additional data and code to make plots using data from The Cancer Genome Atlas and The Human Protein Atlas, which are not included here because of their large size.
}

\references{
  Dick, Jeffrey M. (2016) Proteomic indicators of oxidation and hydration state in colorectal cancer. \emph{PeerJ} \bold{4}, e2238. \doi{10.7717/peerj.2238}

  Dick, Jeffrey M. (2017) Chemical composition and the potential for proteomic transformation in cancer, hypoxia, and hyperosmotic stress. \emph{PeerJ} \bold{5}, e3421. \doi{10.7717/peerj.3421}

  Dick, J. M., Yu, M. and Tan, J. (2020) Distinct trends in chemical composition of proteins from metagenomes in redox and salinity gradients. \emph{bioRxiv}. \doi{10.1101/2020.04.01.020008}

  Dick, J. M. (2020) Water as a reactant in the differential expression of proteins in cancer. \emph{(preprint)}
}

\keyword{package}
