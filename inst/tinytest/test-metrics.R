info <- "nH2O() and Zc() give expected results"
# Get nH2O and Zc for a few proteins the "long way" (using functions in CHNOSZ)
library(CHNOSZ)
basis(c("glutamine", "glutamic acid", "cysteine", "H2O", "O2"))
H2O.ref <- protein.basis(1:6)[, "H2O"] / protein.length(1:6)
O2.ref <- protein.basis(1:6)[, "O2"] / protein.length(1:6)
Zc.ref <- ZC(protein.formula(1:6))

# Get nH2O and Zc using functions in canprot
AAcomp <- thermo()$protein[1:6, ]
H2O.calc <- nH2O(AAcomp, "QEC", terminal_H2O = 1)
O2.calc <- nO2(AAcomp, "QEC")
Zc.calc <- Zc(AAcomp)

# Make the tests
expect_equivalent(H2O.ref, H2O.calc)
expect_equivalent(O2.ref, O2.calc)
expect_equivalent(Zc.ref, Zc.calc)

info <- "nH2O() and Zc() also work on _matrices_ with one row"
# Added on 20220107 for canprot 1.1.2 (functionality used in JMDplots::getmetrics)
library(CHNOSZ)
# Get data frames for proteins with identifiers and amino acid composition
AAcomp6 <- thermo()$protein[1:6, ]
# Extract a numeric _matrix_ with amino acid composition
AAmat6 <- as.matrix(AAcomp6[, 6:25])
# Now get just the first row
AAmat1 <- AAmat6[1, , drop = FALSE]

# Make the tests
expect_equivalent(Zc(AAcomp6), Zc(AAmat6))
expect_equivalent(nH2O(AAcomp6), nH2O(AAmat6))
expect_equivalent(Zc(AAmat6)[1], Zc(AAmat1))
expect_equivalent(nH2O(AAmat6)[1], nH2O(AAmat1))