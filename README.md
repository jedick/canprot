[![CRAN](https://www.r-pkg.org/badges/version/canprot)](https://cran.r-project.org/package=canprot)
[![DOI](https://zenodo.org/badge/64122601.svg)](https://zenodo.org/badge/latestdoi/64122601)

# canprot

Chemical metrics of differentially expressed proteins in cancer and cell
culture proteomics experiments. Data files in the package have amino acid
compositions of proteins obtained from UniProt and >250 published lists of up-
and down-regulated proteins in different cancer types and laboratory
experiments. Functions are provided to calculate chemical metrics including
protein length, grand average of hydropathicity (GRAVY), isoelectric point
(pI), carbon oxidation state, and stoichiometric hydration state; the latter
two are described in [Dick et al.
(2020)](https://doi.org/10.5194/bg-17-6145-2020). The vignettes visualize
differences of chemical metrics between up- and down-regulated proteins and
list literature references for all datasets.

Please use this citation for the package: [Dick (2021)](https://doi.org/10.1002/cso2.1007).

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

This also installs other R packages as dependencies (particularly **xtable**, **knitr** and **rmarkdown**, and their dependencies).

### Building vignettes

The main vignette provides a brief introduction to the package.
There are separate analysis vignettes for each dataset; these are not built when the packagke is installed.

Building the analysis vignettes requires [pandoc](https://pandoc.org) to be available on the system.
To compile one of the analysis vignettes and open it in your browser, use the `mkvig()` function in **canprot**, like this:
```R
library(canprot)
mkvig("3D")
```

## Online vignettes

All vignettes can be viewed at <http://chnosz.net/canprot/doc/index.html>.
