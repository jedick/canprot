\encoding{UTF-8}
\name{rankdiff}
\alias{rankdiff}
\title{Weighted Difference of Sums of Ranks}
\description{
  Calculate rank-sum difference between two groups, weighted by the sizes of the groups.
}

\usage{
  rankdiff(rank1, rank2, n1 = NULL, n2 = NULL, as.fraction=TRUE)
}

\arguments{
  \item{rank1}{numeric, ranks in group 1}
  \item{rank2}{numeric, ranks in group 2}
  \item{n1}{numeric, size of group 1}
  \item{n2}{numeric, size of group 2}
  \item{as.fraction}{logical, calculate the fraction of maximum possible difference?}
}

\details{
In a combined ranking of two groups, the comparison of sum of ranks has an easy interpretation only for groups of equal size.
The weighted rank difference is used to compare groups of unequal size.
The weighting ensures that 1) opposite extreme configurations give weighted rank differences with equal magnitudes, and 2) an evenly distributed (interspersed) ranking of the two groups has a weighted rank difference of zero (Dick, 2016).

If \code{n1} and \code{n2} are not given, \code{rank1} and \code{rank2} are interpreted as vectors holding the ranks for the two groups.
If the sizes of the groups are supplied in \code{n1} and \code{n2}, then the single values or higher-dimensional objects in \code{rank1} and \code{rank2} are interpreted as the non-weighted sums of ranks of the two groups.
}

\examples{
# rankings of H and C in H-H-H-H-C-C-C
rankdiff(1:4, 5:7, as.fraction=FALSE)  # 12
rankdiff(1:4, 5:7)  # 1

# rankings of H and C in C-C-C-H-H-H-H
rankdiff(4:7, 1:3, as.fraction=FALSE)  # -12
rankdiff(4:7, 1:3)  # -1

# rankings of H and C in H-C-H-C-H-C-H
rankdiff(c(1, 3, 5, 7), c(2, 4, 6))  # 0
}

\references{
  Dick, Jeffrey M. (2016) Proteomic indicators of oxidation and hydration state in colorectal cancer. \emph{PeerJ} \bold{4}, e2238. \doi{10.7717/peerj.2238}
}

\concept{Statistical functions}
