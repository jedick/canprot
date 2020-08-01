[![CRAN](https://www.r-pkg.org/badges/version/canprot)](https://cran.r-project.org/package=canprot)
[![DOI](https://zenodo.org/badge/64122601.svg)](https://zenodo.org/badge/latestdoi/64122601)

# canprot

Compositional analysis of differentially expressed proteins in cancer and cell
culture proteomics experiments. The data include lists of up- and
down-regulated proteins in different types of cancer (breast, colorectal,
liver, lung, pancreatic, prostate) and laboratory conditions (hypoxia,
hyperosmotic stress, high glucose, 3D cell culture, and proteins secreted in
hypoxia), together with amino acid compositions computed for protein sequences
obtained from UniProt. Functions are provided to calculate compositional metrics
including protein length, carbon oxidation state, and stoichiometric hydration
state. In addition, phylostrata (evolutionary ages) of protein-coding genes are
compiled using data from [Liebeskind et al. (2016)](https://doi.org/10.1093/gbe/evw113) or
[Trigos et al. (2017)](https://doi.org/10.1073/pnas.1617743114). The vignettes contain
plots of compositional differences, phylostrata (for human proteins), and
references for all datasets.

For more information, see two papers in *PeerJ* ([2016](http://doi.org/10.7717/peerj.2238)
and [2017](http://doi.org/10.7717/peerj.3421)).

## Installation from CRAN

```R
install.packages("canprot")
```

## Installation from GitHub

First install the **remotes** package from CRAN:

```R
install.packages("remotes")
```

Then install **canprot** from GitHub:

```R
remotes::install_github("jedick/canprot")
```

This also installs other R packages as dependencies (particularly **CHNOSZ**, **xtable**, **knitr** and **rmarkdown**, and their dependencies).

### Building vignettes

Building the vignettes requires [pandoc](http://pandoc.org/installing.html) to be available on the system.
See rmarkdown's [Install Pandoc](https://CRAN.R-project.org/package=rmarkdown/vignettes/pandoc.html) vignette for tips on installing pandoc.
With all the dependencies available, the vignettes can be compiled and viewed using the `mkvig()` function in **canprot**, like this:
```R
library(canprot)
mkvig("3D")
```

## Online vignettes

The vignettes can be viewed at <http://chnosz.net/canprot/doc/index.html>.

