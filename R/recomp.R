# recompute a pdat object for a different 'basis' 20191207
recomp <- function(pdat, basis = "MTa") {
  uniprot <- pdat$pcomp$uniprot
  aa <- pdat$pcomp$aa
  pcomp <- protcomp(uniprot, basis, aa)
  pdat$pcomp <- pcomp
  pdat
}
