\encoding{UTF-8}
\name{groupplots}
\alias{groupplots}
\alias{mergedplot}
\title{Plot Potential Diagrams for Groups of Datasets}
\description{
  Plot rank difference of chemical affinities for proteins in various datasets and merge the diagrams.
}

\usage{
  groupplots(group="hypoxia_ZC_down", each100=FALSE, res=50, plot.it=TRUE)
  mergedplot(gpresult, each100=FALSE, res=50)
}

\arguments{
  \item{group}{character, description of datasets to include}
  \item{each100}{logical, rescale rank difference of each dataset individually?}
  \item{res}{numeric, grid resolution for plots}
  \item{plot.it}{logical, make plots?}
  \item{gpresult}{list, value returned by \code{groupplots}}
}

\details{
\code{groupplots} makes weighted rank-difference of affinity (potential) diagrams (see \code{\link{rankplot}}) for each dataset found in the specified \code{group}.
\code{group} consists of three parts joined by an underscore: the type of experiment (\samp{colorectal}, \samp{pancreatic}, \samp{hypoxia}, or \samp{osmotic}; see \code{\link{pdat_}}), the distinguishing compositional variable (\samp{ZC} or \samp{H2O}), and the direction of change of that variable (\samp{up} or \samp{down}).

To identify the datasets in any group, compositional summaries for each dataset are read from pre-calculated tables in \code{extdata/summary}.
Datasets are included for which the absolute mean difference of either \samp{ZC} or \samp{H2O} between up- and down-expressed proteins is greater than 0.01 and the other of \samp{ZC} or \samp{H2O} has \emph{p}-value >= 0.05 and abs(\code{\link{CLES}} - 50) < 10.

\code{groupplots} makes calculations over a large range of \logfO2 and \logaH2O in order to encompass the equipotential lines for most datasets.
This way, the positions of the median and interquartiles of the equipotential lines can be calculated accurately for the \code{mergedplot}, which covers a smaller range of \logfO2 and \logaH2O.
}

\examples{
\dontrun{
gpresult <- groupplots("osmotic_H2O_down", res=25)
mergedplot(gpresult, res=25)

# reproduce Figure 3 of Dick, 2017
ZCgroups <- c("colorectal_ZC_up","pancreatic_ZC_up", "hypoxia_ZC_down")
H2Ogroups <- c("colorectal_H2O_up", "pancreatic_H2O_up", "osmotic_H2O_down")
allgroups <- c(ZCgroups, H2Ogroups)
par(mfrow=c(2, 3))
for(group in allgroups) {
  gpresult <- groupplots(group, plot.it=FALSE)
  mergedplot(gpresult)
  title(main=group)
}}
}

\seealso{
  \code{\link{rankplot}}
}

\concept{Thermodynamics plots}
