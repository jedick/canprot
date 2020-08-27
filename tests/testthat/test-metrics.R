context("metrics")

test_that("H2OAA() and ZCAA() give expected results", {
  # Get nH2O and ZC for a few proteins the "long way" (using functions in CHNOSZ)
  basis(c("glutamine", "cysteine", "acetic acid", "H2O", "O2"))
  H2O.ref <- protein.basis(1:6)[, "H2O"] / protein.length(1:6)
  ZC.ref <- ZC(protein.formula(1:6))

  # Get nH2O and ZC using functions in canprot
  AAcomp <- thermo()$protein[1:6, ]
  H2O.calc <- H2OAA(AAcomp, "QCa")
  ZC.calc <- ZCAA(AAcomp)

  # Make the tests
  expect_equivalent(H2O.ref, H2O.calc)
  expect_equivalent(ZC.ref, ZC.calc)
})

