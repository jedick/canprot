# recompute a pdat object for a different 'basis' 20191207
pdat_recomp <- function(pdat, basis = "rQEC") {
  uniprot <- pdat$pcomp$uniprot
  aa <- pdat$pcomp$aa
  pcomp <- protcomp(uniprot, basis, aa)
  pdat$pcomp <- pcomp
  pdat
}
