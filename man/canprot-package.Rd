\name{canprot-package}
\alias{canprot-package}
\docType{package}
\title{Compositional Analysis of Differentially Expressed Proteins}
\description{
  \pkg{canprot} is a package for analysis of the chemical compositions of proteins from their amino acid compositions.
  The package compiles datasets for differentially expressed proteins in cancer and cell culture conditions from over 250 studies.
}

\section{Overview}{
  This package includes datasets for differential expression of proteins in six cancer types (breast, colorectal, liver, lung, pancreatic, prostate), and four cell culture conditions (hypoxia, hyperosmotic stress, secreted proteins in hypoxia, and 3D compared to 2D growth conditions).
  The hyperosmotic stress data are divided into bacteria, archaea (both high- and low-salt experiments) and eukaryotes; the latter are further divided into salt and glucose experiments.
  Nearly all datasets use UniProt IDs; if not given in the original publications they have been added using the UniProt mapping tool (\url{https://www.uniprot.org/mapping/}).

  The vignettes have plots for each cancer type and cell culture condition and references for all data sources used.
  Because of their size, pre-built vignette HTML files are not included with the package; use \code{\link{mkvig}} to compile and view any of the vignettes.

  The functions in this package were originally derived from code used for the papers of Dick (2016 and 2017).
  Updated data compilations and revised metrics, plots, and vignettes were developed for the papers of Dick et al. (2020) and Dick (2020).
}

\examples{
# List the data files for all studies
# (one study can have more than one dataset)
exprdata <- system.file("extdata/expression", package="canprot")
datafiles <- dir(exprdata, recursive=TRUE)
print(datafiles)
# Show the number of data files for each condition
table(dirname(datafiles))
}

\seealso{
  The \pkg{JMDplots} package on GitHub (\url{https://github.com/jedick/JMDplots}) has vignettes showing analysis of data from The Cancer Genome Atlas and The Human Protein Atlas, which are not included here because of the large size of the data files.
}

\references{
Dick, J. M. (2016) Proteomic indicators of oxidation and hydration state in colorectal cancer. \emph{PeerJ} \bold{4}, e2238. \doi{10.7717/peerj.2238}

Dick, J. M. (2017) Chemical composition and the potential for proteomic transformation in cancer, hypoxia, and hyperosmotic stress. \emph{PeerJ} \bold{5}, e3421. \doi{10.7717/peerj.3421}

Dick, J. M., Yu, M. and Tan, J. (2020) Uncovering chemical signatures of salinity gradients through compositional analysis of protein sequences. \emph{Biogeosciences Discussions}. \doi{10.5194/bg-2020-146}

Dick, J. M. (2020) Water as a reactant in the differential expression of proteins in cancer. \emph{bioRxiv}. \doi{10.1101/2020.04.09.035022}
}

\keyword{package}
