\encoding{UTF-8}
\name{lapply_canprot}
\alias{lapply_canprot}
\title{Parallel Computation}
\description{
  Set up and run parallel computations.
}

\usage{
  lapply_canprot(X, FUN, ..., varlist = NULL, min.length = 10)
}

\arguments{
  \item{X}{vector (atomic or list)}
  \item{FUN}{function or character string naming a function}
  \item{...}{additional arguments to pass to ‘FUN’}
  \item{varlist}{character, names of variables to export}
  \item{min.length}{numeric, minimum length of X for running parallel computations}
}

\details{
If \code{length(X)} is less than \code{min.length}, this function simply calls \code{\link{lapply}} given the \code{X}, \code{FUN} and \code{...} arguments, and returns the result.
Otherwise, a cluster is initiatied by loading \pkg{CHNOSZ} and \pkg{canprot} and setting up the respective data environments (\code{\link[CHNOSZ]{thermo}} and \code{\link{canprot}}).
If \code{varlist} is supplied, the named variables are exported to the global environment of each node in the cluster using \code{\link[parallel]{clusterExport}}.
Then, parallel computations in the cluster are run using \code{\link{parLapply}} given the \code{X}, \code{FUN} and \code{...} arguments.

For examples, see the data vignettes and the help for \code{\link{pdat_}}.
}

\concept{Parallel computation}
