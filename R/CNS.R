# canprot/R/CNS.R
# proteomic differences of elemental abundance (C, N, S) per residue
# 20170124 jmd

CNS <- function(pdat) {
  CNS <- data.frame(pdat$pcomp$residue.formula[, c("C", "N", "S")])
  CNS1 <- CNS[!pdat$up2, ]
  CNS2 <- CNS[pdat$up2, ]
  mean1 <- data.frame(t(colMeans(CNS1)))
  mean2 <- data.frame(t(colMeans(CNS2)))
  colnames(mean1) <- paste0(colnames(CNS), ".mean1")
  colnames(mean2) <- paste0(colnames(CNS), ".mean2")
  diff <- mean2 - mean1
  colnames(diff) <- paste0(colnames(CNS), ".diff")
  p.value <- data.frame(t(sapply(seq_along(CNS), function(icol) stats::wilcox.test(CNS1[, icol], CNS2[, icol])$p.value)))
  colnames(p.value) <- paste0(colnames(CNS), ".p.value")
  CLES <- data.frame(t(sapply(seq_along(CNS), function(icol) 100 * CLES(CNS1[, icol], CNS2[, icol]))))
  colnames(CLES) <- paste0(colnames(CNS), ".CLES")
  cbind(data.frame(dataset = pdat$dataset, description = pdat$description,
                   n1 = nrow(CNS1), n2 = nrow(CNS2), mean1, mean2, diff, p.value, CLES))
}
