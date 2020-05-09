\encoding{UTF-8}
\name{CLES}
\alias{CLES}
\title{Common Language Effect Size}
\description{
  Calculate the common language effect size.
}

\usage{
  CLES(x, y)
}

\arguments{
  \item{x}{numeric, data}
  \item{y}{numeric, data}
}

\details{
The common language statistic is defined for continuous data as \dQuote{the probability that a score sampled at random from one distribution will be greater than a score sampled from some other distribution.} (McGraw and Wong, 1992)

This function calculates the fraction of all possible pairings between elements of \code{x} and \code{y} where the difference (\samp{y} value - \samp{x} value) is positive.

The example uses synthetic data; the actual data used by McGraw and Wong is from NCHS (1987, Tables 13 and 14).
}

\examples{
# Generate synthetic data for heights in inches of 18-24
# year-old males and females
height_male <- rnorm(988, 69.7, 2.8)
height_female <- rnorm(1066, 64.3, 2.6)
# The CLES is approximately 0.92 (McGraw and Wong, 1992)
CLES(height_female, height_male)
}

\references{
  McGraw, Kenneth O. and Wong, S. P. (1992) A common language effect size statistic. \emph{Psychological Bulletin} \bold{11}, 361--365. \doi{10.1037/0033-2909.111.2.361}

  National Center for Health Statistics (1987) \emph{Anthropometric Reference Data and Prevalence of Overweight: United States, 1976-1980.} Data from the National Health Survey, Series 11, No. 238. DHHS Publication (PHS) No. 87-1688. U.S. Government Printing Office, Washington, DC. \url{http://www.cdc.gov/nchs/data/series/sr_11/sr11_238.pdf}

}

\concept{Statistics}
