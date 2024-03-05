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

  Two demos are available by running \code{demo("location")} and \code{demo("redoxin")}.
  The first shows \Zc and pI for human proteins with different subcellular locations.
  The localization data was taken from Thul et al. (2017) and is stored in \file{extdata/protein/TAW+17_Table_S6_Validated.csv}.
  The second shows \Zc vs the midpoint reduction potentials of ferredoxin and thioredoxin in spinach and glutaredoxin and thioredoxin in \emph{E. coli}.
  The sequences were obtained from UniProt and is stored in \file{extdata/fasta/redoxin.fasta}.
  The reduction potential data was taken from Åslund et al. (1997) and Hirasawa et al. (1999) and is stored in \file{extdata/fasta/redoxin.csv}.
  This file also lists UniProt IDs and start and stop positions for the protein chains excluding initiator methionines and signal peptides.
}

\references{

Åslund F, Berndt KD, Holmgren A. 1997. Redox potentials of glutaredoxins and other thiol-disulfide oxidoreductases of the thioredoxin superfamily determined by direct protein-protein redox equilibria. \emph{Journal of Biological Chemistry} \bold{272}(49): 30780–30786. \doi{10.1074/jbc.272.49.30780}

Dick JM, Yu M, Tan J. 2020. Uncovering chemical signatures of salinity gradients through compositional analysis of protein sequences. \emph{Biogeosciences} \bold{17}(23): 6145--6162. \doi{10.5194/bg-17-6145-2020}

Dick JM. 2021. Water as a reactant in the differential expression of proteins in cancer. \emph{Computational and Systems Oncology} \bold{1}(1): e1007. \doi{10.1002/cso2.1007}

Hirasawa M, Schürmann P, Jacquot J-P, Manieri W, Jacquot P, Keryer E, Hartman FC, Knaff DB. 1999. Oxidation-reduction properties of chloroplast thioredoxins, ferredoxin:thioredoxin reductase, and thioredoxin \emph{f}-regulated enzymes. \emph{Biochemistry} \bold{38}(16): 5200–5205. \doi{10.1021/bi982783v}

Thul PJ, Åkesson L, Wiking M, Mahdessian D, Geladaki A, Blal HA, Alm T, Asplund A, Björk L, Breckels LM, et al. 2017. A subcellular map of the human proteome. \emph{Science} \bold{356}(6340): eaal3321. \doi{10.1126/science.aal3321}

}
