\encoding{UTF-8}
\name{human}
\alias{human}
\alias{human_aa}
\alias{human_base}
\alias{human_additional}
\alias{human_extra}
\alias{uniprot_updates}
\title{Amino Acid Compositions of Human Proteins}
\description{
  Data for amino acid compositions of proteins and conversion from old to new UniProt IDs.
}

\details{
The amino acid compositions of human proteins are stored in three files under \code{extdata/protein}.
\itemize{
  \item \code{human_base.rds} contains amino acid compositions of canonical isoforms of manually reviewed proteins in the \href{http://www.uniprot.org/}{UniProt} reference human proteome (computed from sequences in \code{UP000005640_9606.fasta.gz}, dated 2016-04-03).
  \item \code{human_additional.rds} contains amino acid compositions of additional proteins (\ifelse{latex}{\cr}{}\code{UP000005640_9606_additional.fasta.gz}) including isoforms and unreviewed sequences. In version 0.1.5, this file was trimmed to include only those proteins that are used in any of the datasets in the package.
  \item \code{human_extra.csv} contains amino acid compositions of other (\dQuote{extra}) proteins used in a dataset but not listed in one of the files above. These proteins may include obsolete, unreviewed, or newer additions to the UniProt database. Not all sequences here are strictly HUMAN (see the \code{organism} column and the \code{ref} column for the reference keys).
}

On loading the package, the individual data files are read and combined, and the result is assigned to the \code{human_aa} object in the \code{human} environment.

As an aid for processing datasets that list old (obsolete) UniProt IDs, the corresponding new (current) IDs are are stored in \code{uniprot_updates}.
These ID mappings have been manually added as needed for individual datasets, and include proteins from humans as well as other organisms.
\code{\link{check_IDs}} performs the conversion of old to new IDs.
}

\format{
\code{human_aa} is a data frame with 25 columns in the format used for amino acid compositions in \pkg{CHNOSZ} (see \code{\link[CHNOSZ]{thermo}}):

\tabular{lll}{
  \code{protein} \tab character \tab Identification of protein\cr
  \code{organism} \tab character \tab Identification of organism\cr
  \code{ref} \tab character \tab Reference key for source of sequence data\cr
  \code{abbrv} \tab character \tab Abbreviation or other ID for protein (e.g. gene name)\cr
  \code{chains} \tab numeric \tab Number of polypeptide chains in the protein\cr
  \code{Ala}\dots\code{Tyr} \tab numeric \tab Number of each amino acid in the protein
}

The \code{protein} column contains UniProt IDs in the format \code{database|accession-isoform}, where \code{database} is most often \samp{sp} (Swiss-Prot) or \samp{tr} (TrEMBL), and \code{isoform} is an optional suffix indicating the isoform of the protein (particularly in the \code{human_additional} file).
}

\seealso{
Amino acid compositions of non-human proteins are stored in directories \code{archaea}, \code{bacteria}, \code{cow}, \code{dog}, \code{mouse}, \code{rat}, and \code{yeast} under \code{extdata/aa}.
These files can be loaded in \code{\link{protcomp}} via the \code{aa_file} argument, which is used e.g. in \code{\link{pdat_osmotic_bact}}.
}

\examples{
# The number of proteins
nrow(get("human_aa", human))
# The number of old to new ID mappings
nrow(get("uniprot_updates", human))
}

\concept{Amino acid composition}
