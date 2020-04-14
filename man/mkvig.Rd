\encoding{UTF-8}
\name{mkvig}
\alias{mkvig}
\title{Compile and View Vignettes from the Command Line}
\description{
This function compiles the indicated vignette and opens it in the browser.
}

\usage{
  mkvig(vig = "colorectal")
}

\arguments{
  \item{vig}{character, name of vignette without \samp{.Rmd} extension}
}

\details{
Starting with version 0.2, in order to reduce package space and check time, pre-built vignettes are not included in the package.
This function was added to compile the vignettes on demand and view them in a browser.

The available vignettes are listed here:
\itemize{
    \item \emph{Cell culture} -- \samp{hypoxia}, \samp{secreted}, \samp{hyperosmotic}, \samp{glucose}, \samp{a3D}
    \item \emph{Cancer} -- \samp{breast}, \samp{colorectal}, \samp{liver}, \samp{lung}, \samp{pancreatic}, \samp{prostate}
}

Note that pandoc (including pandoc-citeproc), as a system dependency of \pkg{rmarkdown}, must be installed.
See \pkg{rmarkdown}'s \samp{pandoc} vignette for installation tips.
}

\seealso{
The vignettes can also be run using e.g. \code{\link{demo}("glucose")}, and through the interactive help system (\code{\link{help.start}} > Packages > canprot > Code demos).
}

\examples{
\dontrun{
mkvig("colorectal")
}
}
