\encoding{UTF-8}
\name{read.fasta}
\alias{read.fasta}
\alias{count.aa}
\alias{aasum}
\title{Functions for Reading FASTA Files}

\description{
  Read protein amino acid composition or sequences from a file and count numbers of amino acids in given sequences.
}

\usage{
  read.fasta(file, iseq = NULL, ret = "count", lines = NULL, 
    ihead = NULL, start = NULL, stop = NULL, type = "protein", id = NULL)
  count.aa(sequence, start = NULL, stop = NULL, type = "protein")
  aasum(AAcomp, abundance = 1, average = FALSE)
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
  \item{sequence}{character, one or more sequences}
  \item{AAcomp}{data frame, amino acid composition(s) of proteins}
  \item{abundance}{numeric, abundances of proteins}
  \item{average}{logical, return the weighted average of amino acid counts?}
}

\details{
\code{read.fasta} is used to retrieve entries from a FASTA file.
Use \code{iseq} to select the sequences to read (the default is all sequences).

The function returns various data formats depending on the value of \code{ret}:
\itemize{
  \item{count}{data frame of amino acid counts}
  \item{seq}{list of sequences}
  \item{fas}{list of all lines in the FASTA file}
}

When \code{ret} is \samp{count}, the header lines of the file are parsed to obtain protein names that are put into the \code{protein} column in the result.
Furthermore, if a UniProt FASTA header is detected (using the regular expression \code{"\\|......\\|.*_"}), the information there (accession, name, organism) is split into the \code{protein}, \code{abbrv}, and \code{organism} columns of the resulting data frame.
this behavior (which may take a while for large files) can be suppressed by supplying protein names in \code{id}.

To speed up processing, if the line numbers of the header lines were previously determined, they can be supplied in \code{ihead}. 
Optionally, the lines of a previously read file may be supplied in \code{lines} (in this case no file is needed so \code{file} should be set to "").

\code{count.aa} is the underlying function that counts the numbers of each amino acid or nucleic-acid base in one or more sequences.
The matching of letters is case-insensitive.
A message is generated if any character in \code{sequence}, excluding spaces, is not one of the single-letter amino acid or nucleobase abbreviations.
\code{start} and/or \code{stop} can be provided to process a fragment of the sequence.
If only one of \code{start} or \code{stop} is present, the other defaults to 1 (\code{start}) or the length of the respective sequence (\code{stop}).

\code{aasum} sums the amino acid compositions in the input \code{AAcomp} data frame.
It only applies to columns with the three-letter abbreviations of amino acids and to a column named \code{chains} (if present).
The amino acid compositions are multiplied by the indicated \code{abundance} after recycling to the number of proteins.
If \code{average} is TRUE the sum is divided by the number of proteins.
The output has one row and the same number of columns as the input; the value in the non-amino acid columns is taken from the first row of the input.
Note that missing values (NA) of amino acid composition are treated as zero.
}

\value{
\code{read.fasta} returns a list of sequences or lines for \code{ret} equal to \samp{seq} or \samp{fas}, respectively.
Otherwise, (for \code{ret} equal to \samp{count}) a data frame with identifying information and amino acid compositions of proteins; specifically these columns: \code{protein}, \code{organism}, \code{ref}, \code{abbrv}, \code{chains}, \code{Ala}, \code{Cys}, \code{Asp}, \code{Glu}, \code{Phe}, \code{Gly}, \code{His}, \code{Ile}, \code{Lys}, \code{Leu}, \code{Met}, \code{Asn}, \code{Pro}, \code{Gln}, \code{Arg}, \code{Ser}, \code{Thr}, \code{Val}, \code{Trp}, \code{Tyr}.

For \code{aasum}, a one-row data frame.
}

\seealso{
Pass the output of \code{read.fasta} to \code{\link[CHNOSZ]{add.protein}} in the CHNOSZ package to set up thermodynamic calculations for proteins.
}

\examples{
## Reading a protein FASTA file
# The path to the file
file <- system.file("extdata/fasta/KHAB17.fasta", package = "canprot")
# Read the sequences, and print the first one
read.fasta(file, ret = "seq")[[1]]
# Count the amino acids in the sequences
aa <- read.fasta(file)
# Calculate protein length (number of amino acids in each protein)
plength(aa)
# Sum the amino acid compositions
aasum(aa)

# Count amino acids in a sequence
count.aa("GGSGG")
# A message is issued for unrecognized characters
count.aa("AAAXXX")
# Count nucleobases in a sequence
bases <- count.aa("ACCGGGTTT", type = "DNA")
}

\concept{Amino acid composition}
