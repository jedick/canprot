context("metrics")

test_that("H2OAA() and ZCAA() give expected results", {
  # Get nH2O and ZC for a few proteins the "long way" (using functions in CHNOSZ)
  library(CHNOSZ)
  basis(c("glutamine", "glutamic acid", "cysteine", "H2O", "O2"))
  H2O.ref <- protein.basis(1:6)[, "H2O"] / protein.length(1:6)
  O2.ref <- protein.basis(1:6)[, "O2"] / protein.length(1:6)
  ZC.ref <- ZC(protein.formula(1:6))

  # Get nH2O and ZC using functions in canprot
  AAcomp <- thermo()$protein[1:6, ]
  H2O.calc <- H2OAA(AAcomp, "QEC")
  O2.calc <- O2AA(AAcomp, "QEC")
  ZC.calc <- ZCAA(AAcomp)

  # Make the tests
  expect_equivalent(H2O.ref, H2O.calc)
  expect_equivalent(O2.ref, O2.calc)
  expect_equivalent(ZC.ref, ZC.calc)
})

test_that("H2OAA() and ZCAA() also work on _matrices_ with one row", {
  # Added on 20220107 for canprot 1.1.2 (functionality used in JMDplots::getmetrics)
  library(CHNOSZ)
  # Get data frames for proteins with identifiers and amino acid composition
  AAcomp6 <- thermo()$protein[1:6, ]
  # Extract a numeric _matrix_ with amino acid composition
  AAmat6 <- as.matrix(AAcomp6[, 6:25])
  # Now get just the first row
  AAmat1 <- AAmat6[1, , drop = FALSE]

  # Make the tests
  expect_equivalent(ZCAA(AAcomp6), ZCAA(AAmat6))
  expect_equivalent(H2OAA(AAcomp6), H2OAA(AAmat6))
  expect_equivalent(ZCAA(AAmat6)[1], ZCAA(AAmat1))
  expect_equivalent(H2OAA(AAmat6)[1], H2OAA(AAmat1))
})
