\encoding{UTF-8}
\name{qdist}
\alias{qdist}
\title{Quantile Distributions for One Dataset}
\description{
  Make a plot showing quantile distributions for up- and down-regulated proteins.
}

\usage{
  qdist(pdat, vars = c("ZC", "nH2O"), show.steps = FALSE)
}

\arguments{
  \item{pdat}{list, output of a \code{\link{pdat_}} function for a single dataset}
  \item{vars}{character, which variables to plot}
  \item{show.steps}{logical, show the steps using \code{\link{plot.ecdf}}?}
}

\details{

This function makes a quantile distribution plot with lines for both up- and down-regulated proteins.
The variable (\code{var}) can be \samp{ZC}, \samp{H2O}, or both (two plots are made for the latter).
The horizontal axis is the variable and the vertical axis is the quantile point.
A solid black line is drawn for the down-regulated proteins, and a dashed red line for the up-regulated proteins.
The median difference is shown by a gray horizontal line drawn between the distributions at the 0.5 quantile point.

}

\seealso{
\code{\link{diffplot}} to plot median differences for multiple datasets.
}

\examples{
# Plot the data of Jimenez et al., 2010 for colorectal cancer
pdat <- pdat_colorectal("JKMF10")
qdist()
}

\references{
Jimenez, C. R. and Knol, J. C. and Meijer, G. A. and Fijneman, R. J. A. (2010) Proteomics of colorectal cancer: Overview of discovery studies and identification of commonly identified cancer-associated proteins and candidate CRC serum markers. \emph{J. Proteomics} \bold{73}, 1873--1895. \doi{10.1016/j.jprot.2010.06.004}
}

\concept{Plotting functions}
