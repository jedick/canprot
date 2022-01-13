[![CRAN](https://www.r-pkg.org/badges/version/canprot)](https://cran.r-project.org/package=canprot)
[![DOI](https://zenodo.org/badge/64122601.svg)](https://zenodo.org/badge/latestdoi/64122601)

# canprot

Chemical analysis of differentially expressed proteins in cancer and cell
culture proteomics experiments. The data files are derived from >250 published
lists of up- and down-regulated proteins in different cancer types (breast,
colorectal, liver, lung, pancreatic, prostate) and laboratory experiments
(hypoxia, hyperosmotic stress, high glucose, 3D cell culture, and proteins
secreted in hypoxia), together with amino acid compositions for proteins
obtained from UniProt. Functions are provided to calculate chemical metrics
such as protein length, carbon oxidation state, and stoichiometric hydration
state following methods described by [Dick et al.
(2020)](https://doi.org/10.5194/bg-17-6145-2020).  The vignettes have plots of
differences of chemical metrics between up- and down-regulated proteins and
literature references for all datasets.

If you use this package, please cite [Dick (2021)](https://doi.org/10.1002/cso2.1007).

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

Building the vignettes requires [pandoc](https://pandoc.org) to be available on the system.
With all the dependencies available, the vignettes can be compiled and viewed using the `mkvig()` function in **canprot**, like this:
```R
library(canprot)
mkvig("3D")
```

## Online vignettes

The vignettes can be viewed at <http://chnosz.net/canprot/doc/index.html>.

