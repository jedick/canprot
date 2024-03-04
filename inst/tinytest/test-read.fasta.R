# Test added on 20230304
info <- "read.fasta() recognizes type, iseq, start, stop arguments"
# Read protein IDs, sequence start/stop positions, and midpoint potentials
data_file <- system.file("extdata/fasta/redoxin.csv", package = "canprot")
dat <- read.csv(data_file)
# Read header lines
fasta_file <- system.file("extdata/fasta/redoxin.fasta", package = "canprot")
# TODO: change canprot::read.fasta() to read.fasta() after CHNOSZ 2.1.0 is superseded on CRAN 20240304
headers <- canprot::read.fasta(fasta_file, type = "header")
# Locate the sequences in the FASTA file
iseqs <- sapply(dat$ID, grep, x = headers)
expect_equal(iseqs, c(5, 11, 6, 2, 8, 7, 1, 9, 10, 3, 4), check.names = FALSE, info = info)
# Loop over proteins
aalist <- lapply(1:nrow(dat), function(i) {
  # Read the amino acid composition of this protein
  canprot::read.fasta(fasta_file, iseq = iseqs[i], start = dat$Start[i], stop = dat$Stop[i])
})
aa <- do.call(rbind, aalist)
expect_equal(gsub(".*\\|", "", aa$protein), dat$ID, info = info)
expect_equal(plength(aa), dat$Stop - dat$Start + 1, info = info)
Zc_ref <- c(-0.223485, -0.09324, -0.126904, -0.139984, -0.166184,
  -0.159664, -0.048673, -0.235294, -0.05036, -0.22028, -0.265018)
expect_equal(round(Zc(aa), 6), Zc_ref)
