\encoding{UTF-8}
\name{read.fasta}
\alias{read.fasta}
\alias{count.aa}
\alias{aasum}
\title{Functions for reading FASTA files}

\description{
  Read protein amino acid composition or sequences from a file and count numbers of amino acids in given sequences.
}

\usage{
  read.fasta(file, iseq = NULL, type = "count", lines = NULL, 
    ihead = NULL, start = NULL, stop = NULL, molecule = "protein", id = NULL)
  count.aa(sequence, start = NULL, stop = NULL, molecule = "protein")
  aasum(AAcomp, abundance = 1, average = FALSE)
}

\arguments{
  \item{file}{character, path to FASTA file}
  \item{iseq}{numeric, which sequences to read from the file}
  \item{type}{character, type of return value (\samp{count}, \samp{sequence}, \samp{lines}, or \samp{headers})}
  \item{lines}{list of character, supply the lines here instead of reading them from file}
  \item{ihead}{numeric, which lines are headers}
  \item{start}{numeric, position in sequence to start counting}
  \item{stop}{numeric, position in sequence to stop counting}
  \item{molecule}{character, type of molecule (\samp{protein}, \samp{DNA}, or \samp{RNA})}
  \item{id}{character, value to be used for \code{protein} in output table}
  \item{sequence}{character, one or more sequences}
  \item{AAcomp}{data frame, amino acid composition(s) of proteins}
  \item{abundance}{numeric, abundances of proteins}
  \item{average}{logical, return the weighted average of amino acid counts?}
}

\details{
\code{read.fasta} is used to retrieve entries from a FASTA file.
Use \code{iseq} to select the sequences to read (the default is all sequences).

The function returns various data formats depending on the value of \code{type}:
\describe{
  \item{\samp{count}}{data frame of amino acid counts}
  \item{\samp{sequence}}{list of sequences}
  \item{\samp{lines}}{list of lines from the FASTA file (including headers)}
  \item{\samp{headers}}{list of header lines from the FASTA file}
}

When \code{type} is \samp{count}, the header lines of the file are parsed to obtain protein names that are put into the \code{protein} column in the result.
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
The values in these columns are multiplied by the indicated \code{abundance} after recycling to the number of proteins.
The values in these columns are then summed; if \code{average} is TRUE then the sum is divided by the number of proteins.
Proteins with missing values (NA) of amino acid composition or abundance are omitted from the calculation.
The output has one row and the same number of columns as the input; the value in the non-amino acid columns is taken from the first row of the input.
}

\value{
\code{count.aa} returns a data frame with these columns (for proteins): \code{Ala}, \code{Cys}, \code{Asp}, \code{Glu}, \code{Phe}, \code{Gly}, \code{His}, \code{Ile}, \code{Lys}, \code{Leu}, \code{Met}, \code{Asn}, \code{Pro}, \code{Gln}, \code{Arg}, \code{Ser}, \code{Thr}, \code{Val}, \code{Trp}, \code{Tyr}.
For \samp{DNA}, the columns are changed to \code{A}, \code{C}, \code{G}, code{T}, and for \samp{RNA}, the columns are changed to \code{A}, \code{C}, \code{G}, code{U}.

\code{read.fasta} returns a list of sequences (for \code{type} equal to \samp{sequence}) or a list of lines (for \code{type} equal to \samp{lines} or \samp{headers}).
Otherwise, (for \code{type} equal to \samp{count}) a data frame with these columns: \code{protein}, \code{organism}, \code{ref}, \code{abbrv}, \code{chains}, and the columns described above for \code{count.aa}.

\code{aasum} returns a one-row data frame.
}

\seealso{
Pass the output of \code{read.fasta} to \code{\link[CHNOSZ]{add.protein}} in the CHNOSZ package to set up thermodynamic calculations for proteins.
}

\examples{
## Reading a protein FASTA file
# The path to the file
file <- system.file("extdata/fasta/KHAB17.fasta", package = "canprot")
# Read the sequences, and print the first one
read.fasta(file, type = "seq")[[1]]
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
bases <- count.aa("ACCGGGTTT", molecule = "DNA")
}

\concept{Amino acid composition}
