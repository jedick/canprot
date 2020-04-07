\encoding{UTF-8}
\name{canprot}
\alias{canprot}
\alias{human_aa}
\alias{human_base}
\alias{human_additional}
\alias{human_extra}
\alias{uniprot_updates}
\title{Amino Acid Compositions of Human Proteins}
\description{
  Data for amino acid compositions, and updates to UniProt IDs.
}

\details{
The amino acid compositions of human proteins were computed from amino acid sequences in the \href{http://www.uniprot.org/}{UniProt} reference proteome, dated 2016-04-03.
The amino acid compositions of human proteins are stored in three files under \code{extdata/protein}.
\itemize{
  \item \code{human_base.rds} contains amino acid compositions of proteins in the UniProt reference proteome (\code{UP000005640_9606.fasta.gz} containing canonical, manually reviewed sequences).
  \item \code{human_additional.rds} contains amino acid compositions of additional proteins (\ifelse{latex}{\cr}{}\code{UP000005640_9606_additional.fasta.gz}) including isoforms and unreviewed sequences. This file has been trimmed to include only those proteins that are used in the package as of version 0.1.4-40 (2020-01-01).
  \item \code{human_extra.csv} contains amino acid compositions of other (\dQuote{extra}) proteins identified in specific proteomic experiments but not listed in one of the files above; these have been manually added along with the datasets, and may be obsolete, unreviewed, or are newer additions to the UniProt reference proteome. Not all sequences here are strictly HUMAN (see the \code{organism} column), but use the IDs or accession numbers reported in the sources (see the \code{ref} column for the references).
}

On loading the package, the individual data files are read and combined, and the result is assigned to the \code{human_aa} object in the \code{canprot} environment.

As an aid for processing datasets that include old (obsoleted) UniProt IDs, the corresponding new (current) IDs are are stored in \code{uniprot_updates}.
The \code{source} column in \code{uniprot_updates} gives the reference(s) for the datasets that use these IDs.
Conversion of old to new IDs is performed with \code{\link{check_IDs}}.

Amino acid compositions of non-human proteins are stored under \code{extdata/aa} (directories \code{archaea}, \code{bacteria}, \code{yeast}, \code{mouse}, \code{rat}).
Datasets that use these proteins can be identified by the \code{aa_file} argument to \code{\link{protcomp}} in the various \code{\link{pdat_}} functions.
}

\format{
\code{human_aa} follows the format used for storing amino acid compositions in \pkg{CHNOSZ} (see \code{\link[CHNOSZ]{thermo}}):

\tabular{lll}{
  \code{protein} \tab character \tab Identification of protein\cr
  \code{organism} \tab character \tab Identification of organism\cr
  \code{ref} \tab character \tab Reference key for source of compositional data\cr
  \code{abbrv} \tab character \tab Abbreviation or other ID for protein\cr
  \code{chains} \tab numeric \tab Number of polypeptide chains in the protein\cr
  \code{Ala}\dots\code{Tyr} \tab numeric \tab Number of each amino acid in the protein
}

Here, the \code{protein} column contains the UniProt ID (accession), possibly with a suffix indicating the isoform of the protein (esp. from \code{human_additional}).
}

\examples{
nrow(get("human_aa", canprot))
}

\concept{Amino acid composition}
