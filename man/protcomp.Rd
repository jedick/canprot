\encoding{UTF-8}
\name{protcomp}
\alias{protcomp}
\title{Protein Compositions}
\description{
  Get amino acid and chemical compositions of proteins.
}

\usage{
  protcomp(uniprot = NULL, ip = NULL, basis = "QEC",
           aa_file = NULL, updates_file = NULL)
}

\arguments{
  \item{uniprot}{character, UniProt IDs of proteins}
  \item{ip}{numeric, indices of active proteins in CHNOSZ}
  \item{basis}{character, keyword for basis species to use}
  \item{aa_file}{character, file name}
  \item{updates_file}{character, file name}
}

\details{
This function retrieves the amino acid compositions of one or more proteins specified by \code{uniprot} or \code{ip}, then calculates some chemical compositional properties using functions provided by \pkg{CHNOSZ}.
The \code{basis} argument is used to select the basis species using a keyword (see \code{\link[CHNOSZ]{basis}}).
For example, use \samp{CHNOS} for \CO2, \NH3, \H2S, \H2O, and \O2, or \samp{QEC} (the default) for glutamine, glutamic acid, cysteine, \H2O, and \O2.

This function depends on the amino acid compositions of human proteins, which are stored in the \code{\link{canprot}} environment when the package is attached.
If \code{aa_file} is specified, additional amino acid compositions to be considered are read from this file, which should be in the same format as e.g. \code{\link{human_extra}.csv} (see also \code{\link[CHNOSZ]{thermo}$protein}).
If \code{updates_file} is specified, these ID mappings are included alongside the packaged \code{\link{uniprot_updates}}.
}

\value{
The function returns a list with elements \code{protein.formula} (elemental compositions of the proteins), \code{ZC} (average oxidation state of carbon), \code{protein.basis} (compositions of the proteins in terms of the basis species), \code{protein.length} (lengths of the amino acid sequences), \code{residue.basis} (per-residue compositions of the proteins in terms of the basis species), \code{residue.formula} (per-residue elemental compositions of the proteins), and \code{aa} (amino acid compositions of the proteins).
}

\seealso{
\code{\link{cleanup}}
}

\examples{
library(CHNOSZ)
protcomp("P24298")
}

\concept{Chemical composition}
