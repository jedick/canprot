\encoding{UTF-8}
\name{mkvig}
\alias{mkvig}
\title{Compile and View Vignettes from the Command Line}
\description{
Compile the indicated vignette and open it in the browser.
}

\usage{
  mkvig(vig = NULL)
}

\arguments{
  \item{vig}{character, name of a vignette without \samp{.Rmd} extension}
}

\details{
In order to reduce package space and check time, pre-built vignettes are not included in the package.
\code{mkvig} is a convenience function to compile the vignettes on demand and view them in a browser.

The available vignettes for \code{mkvig} are listed here:
\itemize{
    \item \emph{Cell culture} -- \samp{hypoxia}, \samp{secreted}, \samp{osmotic_bact}, \samp{osmotic_euk}, \samp{osmotic_halo}, \samp{glucose}, \samp{3D}
    \item \emph{Cancer} -- \samp{breast}, \samp{colorectal}, \samp{liver}, \samp{lung}, \samp{pancreatic}, \samp{prostate}
}

Note that pandoc (including pandoc-citeproc), as a system dependency of \pkg{rmarkdown}, is required to build the vignettes.
}

\seealso{
The vignettes can also be run using e.g. \code{\link{demo}("glucose")}, and through the interactive help system (\code{\link{help.start}} > Packages > canprot > Code demos).
}

\examples{
\dontrun{
mkvig("colorectal")
}
}

\concept{Vignette utilities}
