\encoding{UTF-8}
\name{get_pdat}
\alias{get_pdat}
\title{Get Protein Data}
\description{
  Get data on protein expression and chemical composition.
}

\usage{
  get_pdat(dataset = NULL, pdat_fun = NULL, basis = "AA")
}

\arguments{
  \item{dataset}{character, specifies which dataset to retrieve}
  \item{pdat_fun}{character, names of additional \code{pdat_} functions}
  \item{basis}{character, keyword for basis species to use}
}

\details{
\code{get_pdat} serves as a wrapper to the various \code{\link{pdat_}} functions.
It queries those functions for the availability of \code{dataset}; if one function provides it, the corresponding dataset is returned.
Zero or more user-defined functions can be named in \code{pdat_fun}.
Leave \code{dataset} as \code{NULL} to list all datasets provided by all \code{pdat_} functions in the package (and by any given in \code{pdat_fun}).
}

\examples{
data(canprot)
get_pdat()
get_pdat("JKMF10")  # same result as pdat_CRC("JKMF10")
}

\concept{Protein expression}
