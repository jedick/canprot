info <- "add.cld() returns letters"
# Modified from ?add.cld
aa <- get("human_aa", canprot)
aa <- aa[!plength(aa) < 20, ]
Zc <- Zc(aa)
ilo <- Zc < -0.15
ihi <- Zc > -0.10
imid <- !ilo & !ihi
nH2O <- nH2O(aa)
nH2Olist <- list(lo.Zc = nH2O[ilo], mid.Zc = nH2O[imid], hi.Zc = nH2O[ihi])
bp <- boxplot(nH2Olist, ylab = cplab$nH2O)
cld <- add.cld(nH2Olist, bp)
expect_equal(cld$letters, c("c", "b", "a"), check.names = FALSE, info = info)
