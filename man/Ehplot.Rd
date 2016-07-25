\encoding{UTF-8}
\name{Ehplot}
\alias{Ehplot}
\title{Plot Eh}
\description{
  Show redox potential (Eh) scale.
}

\usage{
  Ehplot(basis = "AA", T = 37, pH = 7.4)
}

\arguments{
  \item{basis}{character, keyword for basis species to use}
  \item{T}{numeric, temperature in degrees Celsius}
  \item{pH}{numeric, pH}
}

\details{
This function plots selected values of Eh (redox potential) as a function of \logfO2 and \logaH2O.
The lines are labeled with the Eh value in volts.
The temperature and pH can be adjusted using the \code{T} and \code{pH} arguments (currently this only affects the lines, but not the positions of labels).
}

\seealso{
  \code{\link[CHNOSZ]{convert}} (from \pkg{CHNOSZ})
}

\examples{
data(thermo, package="CHNOSZ")
Ehplot()
}

\concept{Thermodynamic description}
