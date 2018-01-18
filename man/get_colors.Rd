\encoding{UTF-8}
\name{get_colors}
\alias{get_colors}
\title{Get Colors}
\description{
  Get colors for rank-difference (potential) diagrams.
}

\usage{
  get_colors(x, max50 = FALSE)
}

\arguments{
  \item{x}{numeric values}
  \item{max50}{logical, use most intense color for all values >= 50?}
}

\details{
\code{get_colors} returns a diverging (blue - light gray - red) color scale.
Blue and red colors are associated with negative and positive values, respectively.
The intensity of the color increases with the magnitude of the value.
For accurate representation, the values should be in a percent scale (i.e. the maximum absolute value is not greater than 100).
By default, a value of +/- 100 corresponds to greatest intensity.
Set \code{max50} to TRUE to compress the scale so that greatest intesity is obtained at values of +/- 50 and higher.
}

\seealso{
  \code{\link{rankplot}}
}

\concept{Protein expression}
