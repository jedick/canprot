# canprot

Datasets are collected here for differentially (up- and down-) expressed
proteins in cancer. All datasets use UniProt IDs, which have been added if not
present in the original publications. Tables of amino acid compositions of
proteins are provided, and the functions are geared toward the exploration of
chemical compositional differences and thermodynamic descriptions using basis
species.

The manual (i.e. help pages) and vignettes showing calculations on datasets for
colorectal cancer can be accessed at <http://chnosz.net/canprot>.

## Installation

First install the *devtools* package from CRAN:

```R
install.packages("devtools")
```

Then install *canprot* from GitHub:

```R
devtools::install_github("jedick/canprot")
```
