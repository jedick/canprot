\name{canprot-package}
\alias{canprot-package}
\docType{package}
\title{Chemical analysis of proteins}
\description{
  Chemical metrics of proteins are important for understanding biomolecular adaptation to environments.
  \pkg{canprot} computes chemical metrics of proteins from their amino acid compositions.
}

\details{
  \itemize{
    \item \code{\link{read.fasta}} -- start here to read amino acid compositions from FASTA files.
    \item \code{\link{metrics}} -- use these functions to calculate chemical metrics from amino acid compositions.
  }

  Chemical metrics include carbon oxidation state (\Zc) stoichiometric oxidation and hydration state (\nO2 and \nH2O) as described in Dick et al. (2020).
  Other variables that can be calculated are protein length, average molecular weight of amino acid residues, grand average of hydropathy (GRAVY), and isoelectric point (pI).

  Differential expression datasets that were in the package up to version 1.1.2, mainly for the paper by Dick (2021), have been moved to \href{https://github.com/jedick/JMDplots}{JMDplots}.

  A demo is available by running \code{demo("canprot")}.
  This shows \Zc and pI for human proteins with different subcellular locations.
  The localization data is taken from Thul et al. (2017) and is stored in \file{extdata/protein/TAW+17_Table_S6_Validated.csv}.
}

\references{
Dick JM, Yu M, Tan J. 2020. Uncovering chemical signatures of salinity gradients through compositional analysis of protein sequences. \emph{Biogeosciences} \bold{17}(23): 6145--6162. \doi{10.5194/bg-17-6145-2020}

Dick JM. 2021. Water as a reactant in the differential expression of proteins in cancer. \emph{Computational and Systems Oncology} \bold{1}(1): e1007. \doi{10.1002/cso2.1007}

Thul PJ, Åkesson L, Wiking M, Mahdessian D, Geladaki A, Blal HA, Alm T, Asplund A, Björk L, Breckels LM, et al. 2017. A subcellular map of the human proteome. \emph{Science} \bold{356}(6340): eaal3321. \doi{10.1126/science.aal3321}
}
