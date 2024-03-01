\encoding{UTF-8}
\name{add.cld}
\alias{add.cld}
\title{Compact Leter Display}
\description{
  Adds compact letter display (significant difference letters) to a boxplot.
}

\usage{
  add.cld(datlist, bp, dx = NULL, dy = NULL)
}

\arguments{
  \item{datlist}{list, list of data values for different groups}
  \item{bp}{list, output of \code{\link{boxplot}}}
  \item{dx}{numeric, offset for letters in the x direction}
  \item{dy}{numeric, offset for letters in the y direction}
}

\details{
This function adds a compact letter display (cld) to an existing boxplot.
It calculates a one-way ANOVA with \code{\link{aov}} followed by Tukey's Honest Significant Differences with \code{\link{TukeyHSD}}, then obtains the cld with \code{\link[multcompView]{multcompLetters4}}.
The letters are added to the plot at the upper right sides of the bars.
Default values for \code{dx} and \code{dy} are computed from the current plot dimensions; these values can be adjusted if needed.
}

\value{
Invisibly returns a list with \code{dx}, \code{dy}, and \code{letters} (letters used for the cld, in the same order as the groups in \code{datlist}).
}

\examples{
# Are there significant differences of nH2O among human proteins with different Zc?
aa <- get("human_aa", canprot)
# Remove extremely short sequences
aa <- aa[!plength(aa) < 20, ]
Zc <- Zc(aa)
ilo <- Zc < -0.15
ihi <- Zc > -0.10
imid <- !ilo & !ihi
nH2O <- nH2O(aa)
nH2Olist <- list(lo.Zc = nH2O[ilo], mid.Zc = nH2O[imid], hi.Zc = nH2O[ihi])
bp <- boxplot(nH2Olist, ylab = cplab$nH2O)
add.cld(nH2Olist, bp)
# Yes, higher Zc is associated with lower nH2O
}

\concept{Statistical functions}
