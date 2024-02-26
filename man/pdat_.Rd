\encoding{UTF-8}
\name{pdat_}
\alias{pdat_}
\title{Get Protein Expression Data}
\description{
  Get differentially expressed proteins and amino acid compositions.
}

\arguments{
  \item{dataset}{character, dataset name}
}

\details{
The \code{pdat_} functions assemble lists of up- and down-regulated proteins and retrieve their amino acid compositions using \code{\link{protcomp}}.
The result can be used with \code{\link{get_comptab}} to make a table of chemical metrics that can then be plotted with \code{\link{diffplot}}.

If \code{dataset} is \samp{2020} (the default) or \samp{2017}, the function returns the names of all datasets in the compilation for the respective year.

Each dataset name starts with a reference key indicating the study (i.e. paper or other publication) where the data were reported.
The reference keys are made by combining the first characters of the authors' family names with the 2-digit year of publication.

If a study has more than one dataset, the reference key is followed by an underscore and an identifier for the particular dataset.
This identifier is saved in the variable named \code{stage} in the functions, but can be any descriptive text.

To retrieve the data, provide a single dataset name in the \code{dataset} argument.
Protein expression data is read from the CSV files stored in \code{extdata/expression/}, under the subdirectory corresponding to the name of the \code{pdat_} function.
Some of the functions also read amino acid compositions (e.g. for non-human proteins) from the files in \code{extdata/aa/}.

Descriptions for each function:

}

\value{
A list consisting of:
\describe{
  \item{\code{dataset}}{Name of the dataset}
  \item{\code{description}}{Descriptive text for the dataset, used for making the tables in the vignettes}
  \item{\code{pcomp}}{UniProt IDs together with amino acid compositions obtained using \code{\link{protcomp}}}
  \item{\code{up2}}{Logical vector with length equal to the number of proteins; TRUE for up-regulated proteins and FALSE for down-regulated proteins}
}
}

\examples{
## List datasets in the 2017 complilation for colorectal cancer
#pdat_colorectal(2017)
## Get proteins and amino acid compositions for one dataset
#pdat_colorectal("JKMF10")
}

\concept{Protein data}
