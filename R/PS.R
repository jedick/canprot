# canprot/PS.R
# get phylostrata for given UniProt IDs
# 20191127 jmd

PS <- function(uniprot) {
  file <- system.file("extdata/phylostrata/TPPG17.csv", package = "canprot")
  dat <- read.csv(file, as.is = TRUE)
  iPS <- match(uniprot, dat$Entry)
  dat$Phylostrata[iPS]
}
