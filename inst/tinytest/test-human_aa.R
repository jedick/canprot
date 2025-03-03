# Tests added on 20250303

info <- "Retrieval using UniProt ID"
pdat <- human_aa("P24298")
expect_equal(pdat$protein, "sp|P24298|ALAT1", info = info)

info <- "Warning for duplicates"
expect_warning(human_aa(rep("P24298", 2), warn_if_duplicated = TRUE), "some uniprot IDs are duplicated: P24298", info = info)

info <- "Get amino acid compositions from file"
# NOTE: this file doesn't contain human proteins, but we can use it anyway
aa_file <- system.file("extdata/aa/methanogen_aa.csv", package = "canprot")
pdat <- human_aa("2190", aa_file = aa_file)
expect_equal(pdat$organism, "Methanocaldococcus jannaschii")
