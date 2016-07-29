# canprot

This is a package for [R](http://r-project.org).

Datasets are collected here for differentially (up- and down-) expressed
proteins in cancer. All datasets use UniProt IDs, which have been added if not
present in the original publications. Tables of amino acid compositions of
proteins are provided, and the functions are geared toward the exploration of
chemical compositional differences and thermodynamic descriptions using basis
species.

The manual (i.e. help pages) and vignettes, which show calculations using
datasets for colorectal cancer, can be accessed at
<http://chnosz.net/canprot/html/00Index.html>.

## Installation

First install the **devtools** package from CRAN:

```R
install.packages("devtools")
```

Then install **canprot** from GitHub:

```R
devtools::install_github("jedick/canprot")
```

## Installation with vignettes

To install the package including the vignettes:

```R
devtools::install_github("jedick/canprot", build_vignettes=TRUE)
```

You may have to re-run this command one or more times. Note that this pulls in
more R packages as dependencies, and you also have to have
[pandoc](http://pandoc.org/installing.html) installed.
