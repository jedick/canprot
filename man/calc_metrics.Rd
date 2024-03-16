\encoding{UTF-8}
\name{calc_metrics}
\alias{calc_metrics}
\title{Calculate one or more chemical metrics}
\description{
Calculates selected chemical metrics from amino acid composition(s) of proteins.
}

\usage{
  calc_metrics(AAcomp, metrics = c("Zc", "nO2", "nH2O"), ...)
}

\arguments{
  \item{AAcomp}{data frame with amino acid compositions}
  \item{metrics}{character, chemical metrics to calculate}
  \item{...}{additional arguments passed to individual functions}
}

\details{

This is a wrapper function for the functions described at \code{\link{metrics}}.
This wrapper implements case-insensitve matching to the function names (e.g., \code{zc} matches \code{\link{Zc}}).
The following additional shortcuts are defined (e.g., \code{length} can be used in place of the actual function name, \code{\link{plength}}):

\tabular{ll}{
  \code{length} \tab \code{\link{plength}} \cr
  \code{H/C}, \code{H_C} \tab \code{\link{HC}} \cr
  \code{N/C}, \code{N_C} \tab \code{\link{NC}} \cr
  \code{O/C}, \code{O_C} \tab \code{\link{OC}} \cr
  \code{S/C}, \code{S_C} \tab \code{\link{SC}} \cr
}

}

\seealso{
\code{\link{metrics}}
}

\value{
A data frame with the same number of rows as \code{AAcomp} and one column of numeric values for each of the \code{metrics}.
An error is produced if any of the \code{metrics} is not available for calculation.
}

\examples{
# Define the acid composition of alanylglycine
AG <- data.frame(Ala = 1, Gly = 1)
# Calculate default metrics (Zc, nO2, nH2O)
calc_metrics(AG)
# Calculate selected metrics
calc_metrics(AG, c("H/C", "O/C", "Length"))
}
