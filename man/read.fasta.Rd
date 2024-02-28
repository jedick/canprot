\encoding{UTF-8}
\name{read.fasta}
\alias{read.fasta}
\alias{count.aa}
\title{Functions for Reading FASTA Files}

\description{
  Read protein amino acid composition or sequences from a file and count numbers of amino acids in given sequences.
}

\usage{
  read.fasta(file, iseq = NULL, ret = "count", lines = NULL, 
    ihead = NULL, start = NULL, stop = NULL, type = "protein", id = NULL)
  count.aa(seq, start = NULL, stop = NULL, type = "protein")
}

\arguments{
  \item{file}{character, path to FASTA file}
  \item{iseq}{numeric, which sequences to read from the file}
  \item{ret}{character, specification for type of return (count, sequence, or FASTA format)}
  \item{lines}{list of character, supply the lines here instead of reading them from file}
  \item{ihead}{numeric, which lines are headers}
  \item{start}{numeric, position in sequence to start counting}
  \item{stop}{numeric, position in sequence to stop counting}
  \item{type}{character, sequence type (protein or DNA)}
  \item{id}{character, value to be used for \code{protein} in output table}
  \item{seq}{character, amino acid sequence of a protein}
}

\details{
\code{read.fasta} is used to retrieve entries from a FASTA file.
Use \code{iseq} to select the sequences to read (the default is all sequences).
The function returns various formats depending on the value of \code{ret}.
The default \samp{count} returns a data frame of amino acid counts, \samp{seq} returns a list of sequences, and \samp{fas} returns a list of lines extracted from the FASTA file, including the headers (this can be used e.g. to generate a new FASTA file with only the selected sequences).
If the line numbers of the header lines were previously determined, they can be supplied in \code{ihead}. 
Optionally, the lines of a previously read file may be supplied in \code{lines} (in this case no file is needed so \code{file} should be set to "").
When \code{ret} is \samp{count}, the names of the proteins in the resulting data frame are parsed from the header lines of the file, unless \code{id} is provided.
If \code{id} is not given, and a UniProt FASTA header is detected (regular expression \code{"\\|......\\|.*_"}), information there (accession, name, organism) is split into the \code{protein}, \code{abbrv}, and \code{organism} columns of the resulting data frame.

\code{count.aa} counts the occurrences of each amino acid or nucleic-acid base in a sequence (\code{seq}).
The matching of letters is case-insensitive.
A warning is generated if any character in \code{seq}, excluding spaces, is not one of the single-letter amino acid or nucleobase abbreviations.
\code{start} and/or \code{stop} can be provided to count a fragment of the sequence (extracted using \code{\link{substr}}).
If only one of \code{start} or \code{stop} is present, the other defaults to 1 (\code{start}) or the length of the sequence (\code{stop}).
}

\value{
\code{read.fasta} returns a list of sequences or lines for \code{ret} equal to \samp{seq} or \samp{fas}, respectively.
Otherwise, (for \code{ret} equal to \samp{count}) a data frame with identifying information and amino acid compositions of proteins; specifically these columns: \code{protein}, \code{organism}, \code{ref}, \code{abbrv}, \code{chains}, \code{Ala}, \code{Cys}, \code{Asp}, \code{Glu}, \code{Phe}, \code{Gly}, \code{His}, \code{Ile}, \code{Lys}, \code{Leu}, \code{Met}, \code{Asn}, \code{Pro}, \code{Gln}, \code{Arg}, \code{Ser}, \code{Thr}, \code{Val}, \code{Trp}, \code{Tyr}.
}

\seealso{
Pass the output of \code{read.fasta} to \code{\link[CHNOSZ]{add.protein}} in the CHNOSZ package to set up thermodynamic calculations for proteins.
}

\examples{
## Reading a protein FASTA file
# The path to the file
file <- system.file("extdata/fasta/EF-Tu.aln", package = "canprot")
# Read the sequences, and print the first one
read.fasta(file, ret = "seq")[[1]]
# Count the amino acids in the sequences
aa <- read.fasta(file)
# Calculate protein length (number of amino acids in each protein)
plength(aa)

# Count amino acids in a sequence
count.aa("GGSGG")
# A message is issued for unrecognized characters
count.aa("AAAXXX")
# Count nucleobases in a sequence
bases <- count.aa("ACCGGGTTT", type = "DNA")
}
