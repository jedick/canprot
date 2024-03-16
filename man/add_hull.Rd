\encoding{UTF-8}
\name{add_hull}
\alias{add_hull}
\title{Add convex hull}
\description{
  Adds a convex hull around the data.
}

\usage{
  add_hull(x, y = NULL, ...)
}

\arguments{
  \item{x}{x values}
  \item{y}{y values}
  \item{...}{arguments for \code{\link{polygon}}}
}

\details{
\code{add_hull} draws a convex hull around the points given in \code{x}, \code{y}.
This function is a wrapper for \code{\link{chull}} and \code{\link{polygon}}.
}

\value{
Invisibly returns the result from \code{chull}.
}

\examples{
dat <- iris[, 1:2]
plot(dat)
add_hull(dat)
}

\concept{Graphics functions}
