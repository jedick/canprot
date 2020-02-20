[![DOI](https://zenodo.org/badge/64122601.svg)](https://zenodo.org/badge/latestdoi/64122601)

# canprot

Compositional analysis of differentially expressed proteins in cancer and cell
culture proteomics experiments. The data include lists of up- and
down-regulated proteins in different types of cancer (colorectal, pancreatic,
breast, lung, prostate) and laboratory conditions (hypoxia, hyperosmotic
stress, 3D cell culture, and proteins secreted in hypoxia), together with amino
acid compositions calculated for protein sequences from UniProt. Functions are
provided to compute compositional metrics including protein length, carbon
oxidation state, and stoichiometric hydration state. In addition, phylostrata
(evolutionary ages) of protein-coding genes are available using data from
[Liebeskind et al. (2016)](https://doi.org/10.1093/gbe/evw113) or [Trigos et
al. (2017)](https://doi.org/10.1073/pnas.1617743114). The vignettes show
references and plots of differences of chemical composition and phylostrata for
all datasets.

For more information, see two papers in *PeerJ* ([2016](http://doi.org/10.7717/peerj.2238)
and [2017](http://doi.org/10.7717/peerj.3421)).

The manual (help pages) and vignettes can be viewed at
<http://chnosz.net/canprot/html/00Index.html>.

## Installation from CRAN

```R
install.packages("canprot")
```

## Installation from Github

First install the **remotes** package from CRAN:

```R
install.packages("remotes")
```

Then install **canprot** from Github:

```R
remotes::install_github("jedick/canprot")
```

This also installs the **CHNOSZ** and **xtable** packages.

### Building vignettes

Building the vignettes requires [pandoc](http://pandoc.org/installing.html) on the system.
To install the package with the vignettes:

```R
remotes::install_github("jedick/canprot", dependencies = TRUE, build_vignettes = TRUE)
```

This installs more R packages as dependencies (particularly **knitr** and **rmarkdown**).
