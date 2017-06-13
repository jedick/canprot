\encoding{UTF-8}
\name{get_pdat}
\alias{get_pdat}
\title{Get Protein Data}
\description{
  Get data on protein expression and chemical composition.
}

\usage{
  get_pdat(dataset = NULL, pdat_fun = "pdat_CRC", basis = "QEC")
}

\arguments{
  \item{dataset}{character, specifies which dataset to retrieve}
  \item{pdat_fun}{character, names of additional \code{pdat_} functions}
  \item{basis}{character, keyword for basis species to use}
}

\details{
\code{get_pdat} serves as a wrapper to the various \code{\link{pdat_}} functions.
The names of the functions, which can include any user-defined functions in addition to those in the package, are specified in \code{pdat_fun}.
Those functions are queried for the availability of \code{dataset}; if one function provides it, the corresponding dataset is returned.
Leave \code{dataset} as \code{NULL} to list all datasets provided by the \code{pdat_} functions.
}

\examples{
data(canprot)
get_pdat()
get_pdat("JKMF10")  # same result as pdat_CRC("JKMF10")
}

\concept{Protein expression}
