[![DOI](https://zenodo.org/badge/64122601.svg)](https://zenodo.org/badge/latestdoi/64122601)

# canprot

Datasets are collected here for differentially (up- and down-)
expressed proteins identified in proteomic studies of cancer and in cell
culture experiments. Tables of amino acid compositions of proteins are
used for calculations of chemical composition, projected into selected
basis species. Plotting functions are used to visualize the compositional
differences and thermodynamic potentials for proteomic transformations.

For more information, see two papers in *PeerJ* ([2016](http://doi.org/10.7717/peerj.2238)
and [2017](http://doi.org/10.7717/peerj.3421)).

The manual (help pages) and vignettes can be viewed at
<http://chnosz.net/canprot/html/00Index.html>.

## Installation from CRAN

```R
install.packages("canprot")
```

## Installation from Github

First install the **devtools** package from CRAN:

```R
install.packages("devtools")
```

Then install **canprot** from Github:

```R
devtools::install_github("jedick/canprot")
```

## Building vignettes

To install the package including the vignettes:

```R
devtools::install_github("jedick/canprot", build_vignettes = TRUE)
```

You may need to re-run this command one or more times. Note that this pulls in
more R packages as dependencies, and [pandoc](http://pandoc.org/installing.html)
is also required.
