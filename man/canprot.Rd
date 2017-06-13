\encoding{UTF-8}
\name{canprot}
\docType{data}
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
These amino acid compositions were compiled from amino acid sequences downloaded from \href{http://www.uniprot.org/}{UniProt}.
Amino acid sequences of human proteins were obtained from files in the UniProt reference proteome, dated 2016-04-03,
downloaded from \url{ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/}. 

The amino acid compositions of human proteins are stored in three files.
\code{human_base.Rdata} contains amino acid compositions of proteins in the UniProt reference proteome (\code{UP000005640_9606.fasta.gz} containing canonical, manually reviewed sequences).
\code{human_additional.Rdata} contains amino acid compositions of additional proteins in the UniProt reference proteome (\ifelse{latex}{\cr}{}\code{UP000005640_9606_additional.fasta.gz} containing isoforms and unreviewed sequences).
\code{human_extra.csv} contains amino acid compositions of other (\dQuote{extra}) proteins identified in proteomic experiments but not listed in one of the files above.

On loading the data (with \code{data(canprot)}), the individual data files are read and combined using \code{\link{rbind}}, and the result is assigned to the \code{human_aa} object in the \code{canprot} environment.

As an aid for processing some datasets that use old (obsoleted) UniProt IDs, the corresponding new (current) IDs are are stored in \code{uniprot_updates}.
\code{uniprot_updates} also lists the source (i.e. reference key) that uses each old ID.
}

\usage{data(canprot)}

\format{
The columns of \code{human_aa} are compatible with the layout used for amino acid compositions in \pkg{CHNOSZ} (see \code{\link[CHNOSZ]{thermo}}):

\tabular{lll}{
  \code{protein} \tab character \tab Identification of protein\cr
  \code{organism} \tab character \tab Identification of organism\cr
  \code{ref} \tab character \tab Reference key for source of compositional data\cr
  \code{abbrv} \tab character \tab Abbreviation or other ID for protein\cr
  \code{chains} \tab numeric \tab Number of polypeptide chains in the protein\cr
  \code{Ala}\dots\code{Tyr} \tab numeric \tab Number of each amino acid in the protein
}

Here, the \code{protein} column contains the UniProt ID (accession), possibly with a suffix indicating the isoform of the protein (esp. from \code{human_additional.csv}).
}

\examples{
data(canprot)
nrow(human_aa)
}

\concept{Amino acid composition}
