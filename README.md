<!-- badges: start -->
[![CRAN](https://img.shields.io/badge/dynamic/yaml?url=https%3A%2F%2Fcloud.r-project.org%2Fweb%2Fpackages%2Fcanprot%2FDESCRIPTION&query=%24.Version&logo=r&label=CRAN&color=4bc51e)](https://cran.r-project.org/package=canprot)
[![DOI](https://zenodo.org/badge/64122601.svg)](https://zenodo.org/badge/latestdoi/64122601)
<!-- badges: end -->

# canprot

Chemical analysis of proteins based on their amino acid compositions. Amino
acid compositions can be read from FASTA files and used to calculate chemical
metrics including carbon oxidation state and stoichiometric hydration state as
described in [Dick et al.  (2020)](https://doi.org/10.5194/bg-17-6145-2020).
Other properties that can be calculated are protein length, grand average of
hydropathicity (GRAVY), isoelectric point (pI), and molecular weight (MW). A
database of amino acid compositions of human proteins derived from UniProt is
provided.

Please use this citation for the package: [Dick (2021)](https://doi.org/10.1002/cso2.1007).

## Installation of developemnt version from GitHub

This version has the `read.fasta()` function that was previously in [CHNOSZ](https://github.com/jedick/CHNOSZ).

First install the **remotes** package from CRAN:

```R
install.packages("remotes")
```

Then install **canprot** from GitHub.
This also installs several other R packages as dependencies:

```R
remotes::install_github("jedick/canprot")
```

## Installation of released version from CRAN

This version does not have all the features in the development version.

```R
install.packages("canprot")
```

### Guide to vignettes

All vignettes can be viewed at <https://chnosz.net/canprot/vignettes>.
The main vignette provides a brief introduction to the package.
