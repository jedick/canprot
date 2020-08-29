\encoding{UTF-8}
\name{Ehplot}
\alias{Ehplot}
\title{Plot Eh}
\description{
  Show redox potential (Eh) scale.
}

\usage{
  Ehplot(T = 37, pH = 7.4,
    xlim = c(-75, -55), ylim = c(-10, 10), dy = 0.2)
}

\arguments{
  \item{T}{numeric, temperature in degrees Celsius}
  \item{pH}{numeric, pH}
  \item{xlim}{numeric, x axis limits}
  \item{ylim}{numeric, y axis limits}
  \item{dy}{numeric, y spacing between iso-\logfO2 lines}
}

\details{
This function plots selected values of Eh (redox potential) as a function of \logfO2 and \logaH2O.
The lines are labeled with the Eh value in volts.
The temperature and pH can be adjusted using the \code{T} and \code{pH} arguments; this affects the lines, but not the positions of labels.
}

\seealso{
This is used in \code{inst/oldvignettes/colorectal.Rmd} to reproduce Figure 6 of Dick (2016).
Other conversions involving Eh are possible with \code{\link[CHNOSZ]{convert}} (from \pkg{CHNOSZ}).
}

\examples{
Ehplot()
}

\references{
  Dick, J. M. (2016) Proteomic indicators of oxidation and hydration state in colorectal cancer. \emph{PeerJ} \bold{4}, e2238. \doi{10.7717/peerj.2238}
}

\concept{Thermodynamics plots}
