\encoding{UTF-8}
\name{qdiff}
\alias{qdiff}
\title{Plot Compositional Differences for One Dataset}
\description{
  Make a plot showing quantiles of compositional metrics for up- and down-regulated proteins.
}

\usage{
  qdiff(pd = pdat_colorectal("JKMF10"), vars = c("ZC", "nH2O"))
}

\arguments{
  \item{pd}{list, output of a \code{\link{pdat_}} function for a single dataset}
  \item{vars}{character, which variables to plot}
}

\details{

This function makes a quantile plot with lines for both up- and down-regulated proteins.
The variable (\code{var}) can be \samp{ZC}, \samp{H2O}, or both (two panels are made).
The horizontal axis is the variable and the vertical axis is the quantile.
The solid black line is for the down-regulated proteins, and the dashed red line is for the up-regulated proteins.
A horizontal line is added to show the median difference (0.5 quantile).

}

\seealso{
\code{\link{qdiff}} to plot median differences for many datasets.
}

\examples{
# uses the data of Jimenez et al., 2020 for colorectal cancer
qdiff()
}

\references{
Jimenez, C. R. and Knol, J. C. and Meijer, G. A. and Fijneman, R. J. A. (2010) Proteomics of colorectal cancer: Overview of discovery studies and identification of commonly identified cancer-associated proteins and candidate CRC serum markers. \emph{J. Proteomics} \bold{73}, 1873--1895. \doi{10.1016/j.jprot.2010.06.004}
}

\concept{Chemical composition}
