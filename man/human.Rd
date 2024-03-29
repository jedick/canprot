\encoding{UTF-8}
\name{human}
\alias{canprot}
\alias{human}
\alias{human.aa}
\alias{human.base}
\alias{human.additional}
\alias{human.extra}
\title{Amino acid compositions of human proteins}
\description{
  Amino acid compositions of human proteins derived from UniProt.
}

\details{

The amino acid compositions of human proteins are stored in three files under \code{extdata/protein}.
\itemize{
  \item \code{human.base.rds} contains amino acid compositions of canonical isoforms of manually reviewed proteins in the \href{https://www.uniprot.org/}{UniProt} reference human proteome (computed from sequences in \code{UP000005640_9606.fasta.gz}, dated 2016-04-03).
  \item \code{human.additional.rds} contains amino acid compositions of additional proteins (\ifelse{latex}{\cr}{}\code{UP000005640_9606_additional.fasta.gz}) including isoforms and unreviewed sequences. In version 0.1.5, this file was trimmed to include only those proteins that are used in any of the datasets in the package.
  \item \code{human.extra.csv} contains amino acid compositions of other (\dQuote{extra}) proteins used in a dataset but not listed in one of the files above. These proteins may include obsolete, unreviewed, or newer additions to the UniProt database. Most, but not all, sequences here are HUMAN (see the \code{organism} column and the \code{ref} column for the reference keys).
}

On loading the package, the individual data files are read and combined, and the result is assigned to the \code{human.aa} object in the \code{canprot} environment.

}

\format{
\code{human.aa} is a data frame with 25 columns in the format used for amino acid compositions in \pkg{CHNOSZ} (see \code{\link[CHNOSZ]{thermo}}):

\tabular{lll}{
  \code{protein} \tab character \tab Identification of protein\cr
  \code{organism} \tab character \tab Identification of organism\cr
  \code{ref} \tab character \tab Reference key for source of sequence data\cr
  \code{abbrv} \tab character \tab Abbreviation or other ID for protein (e.g. gene name)\cr
  \code{chains} \tab numeric \tab Number of polypeptide chains in the protein\cr
  \code{Ala}\dots\code{Tyr} \tab numeric \tab Number of each amino acid in the protein
}

The \code{protein} column contains UniProt IDs in the format \code{database|accession-isoform}, where \code{database} is most often \samp{sp} (Swiss-Prot) or \samp{tr} (TrEMBL), and \code{isoform} is an optional suffix indicating the isoform of the protein (particularly in the \code{human.additional} file).
}

\seealso{
  \code{\link{human_aa}} gets amino acid compositions for human proteins specified by their UniProt IDs.
}

\examples{
# The number of proteins
nrow(get("human.aa", canprot))
}

\concept{Amino acid composition}
