\encoding{UTF-8}
\name{calc.metrics}
\alias{calc.metrics}
\title{Calculate one or more chemical metrics}
\description{
Calculates selected chemical metrics from amino acid composition(s) of one or more proteins.
}

\usage{
  calc.metrics(AAcomp, metrics = c("Zc", "nO2", "nH2O"))
}

\arguments{
  \item{AAcomp}{data frame of amino acid composition with column names being the 3-letter abbreviations of the amino acids (first letter uppercase)}
  \item{metrics}{character, chemical metrics to calculate}
}

\details{

This is a wrapper function for the functions described at \code{\link{metrics}}.
This wrapper implements case-insensitve matching to the function name (so \code{zc} matches \code{\link{Zc}}) and some shortcuts (so \code{length} matches \code{\link{plength}}).

The following metrics are available:
\tabular{ll}{
  \code{Zc} \tab Carbon oxidation state \cr
  \code{nO2} \tab Stoichiometric oxidation state \cr
  \code{nH2O} \tab Stoichiometric hydration state \cr
  \code{GRAVY} \tab Grand average of hydropathy \cr
  \code{pI} \tab Isoelectric point \cr
  \code{MW} \tab Molecular weight \cr
  \code{plength} or \code{length} \tab Protein length \cr
  \code{H/C}, \code{H_C}, or \code{HC} \tab H/C ratio \cr
  \code{N/C}, \code{N_C}, or \code{NC} \tab N/C ratio \cr
  \code{O/C}, \code{O_C}, or \code{OC} \tab O/C ratio \cr
  \code{S/C}, \code{S_C}, or \code{SC} \tab S/C ratio \cr
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
# Calculate default metrics (Zc, nO2, nO2)
calc.metrics(AG)
# Calculate selected metrics
calc.metrics(AG, c("GRAVY", "pI"))
}
