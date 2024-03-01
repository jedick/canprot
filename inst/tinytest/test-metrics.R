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

# Following tests moved here from metrics.Rd 20240229

info <- "Check with Zc of alanine and glycince calculated in CHNOSZ"
Zc.Gly <- CHNOSZ::ZC("C2H5NO2")
Zc.Ala <- CHNOSZ::ZC("C3H7NO2")
# Define the composition of a Gly-Ala-Gly tripeptide
AAcomp <- data.frame(Gly = 2, Ala = 1)
# Calculate the Zc of the tripeptide (value: 0.571)
Zc.GAG <- Zc(AAcomp)
# This is equal to the carbon-number-weighted average of the amino acids
nC.Gly <- 2 * 2
nC.Ala <- 1 * 3
Zc.average <- (nC.Gly * Zc.Gly + nC.Ala * Zc.Ala) / (nC.Ala + nC.Gly)
expect_equal(Zc.GAG, Zc.average)

# Compute the per-residue nH2O of Gly-Ala-Gly
basis("QEC")
nH2O.GAG <- CHNOSZ::species("Gly-Ala-Gly")$H2O
# Divide by the length to get residue average (we keep the terminal H-OH)
nH2O.residue <- nH2O.GAG / 3
# Compare with the value calculated by nH2O() (-0.2)
nH2O.canprot <- nH2O(AAcomp, "QEC", terminal_H2O = 1)
expect_equal(nH2O.residue, nH2O.canprot)

# Added 20240301
info <- "Volume calculation gives expected resuls"
# In CHNOSZ: aa <- pinfo(pinfo("LYSC_CHICK"))
aa <- structure(list(protein = "LYSC", organism = "CHICK", ref = "UniProt", 
    abbrv = "P00698", chains = 1L, Ala = 12, Cys = 8, Asp = 7, 
    Glu = 2, Phe = 3, Gly = 12, His = 1, Ile = 6, Lys = 6, Leu = 8, 
    Met = 2, Asn = 14, Pro = 2, Gln = 3, Arg = 11, Ser = 10, 
    Thr = 7, Val = 6, Trp = 6, Tyr = 3), row.names = 6L, class = "data.frame")
# In CHNOSZ: info(info("LYSC_CHICK"))$V / 129 == 80.78211
expect_equal(as.numeric(V0(aa, terminal_H2O = 1)), 80.78211, tolerance = 1e-4, scale = 1)
