# Test added on 20230304
info <- "read_fasta() recognizes type, iseq, start, stop arguments"
# Read protein IDs, sequence start/stop positions, and midpoint potentials
data_file <- system.file("extdata/fasta/redoxin.csv", package = "canprot")
dat <- read.csv(data_file)
# Read header lines
fasta_file <- system.file("extdata/fasta/redoxin.fasta", package = "canprot")
headers <- read_fasta(fasta_file, type = "header")
# Locate the sequences in the FASTA file
iseqs <- sapply(dat$ID, grep, x = headers)
expect_equal(iseqs, c(5, 11, 6, 2, 8, 7, 1, 9, 10, 3, 4), check.names = FALSE, info = info)
# Loop over proteins
aalist <- lapply(1:nrow(dat), function(i) {
  # Read the amino acid composition of this protein
  read_fasta(fasta_file, iseq = iseqs[i], start = dat$start[i], stop = dat$stop[i])
})
aa <- do.call(rbind, aalist)
expect_equal(gsub(".*\\|", "", aa$protein), dat$ID, info = info)
expect_equal(plength(aa), dat$stop - dat$start + 1, info = info)
Zc_ref <- c(-0.223485, -0.09324, -0.126904, -0.139984, -0.166184,
  -0.159664, -0.048673, -0.235294, -0.05036, -0.22028, -0.265018)
expect_equal(round(Zc(aa), 6), Zc_ref, info = info)

# Test for protein added on 20230308
info <- "read_fasta() handles 0-length 'iseq' argument for proteins"
fasta_file <- system.file("extdata/fasta/redoxin.fasta", package = "canprot")
no_protein_dat <- read_fasta(fasta_file, iseq = numeric(), molecule = "protein")
expect_equal(nrow(no_protein_dat), 0, info = info)
expect_equal(colnames(no_protein_dat)[1], "protein", info = info)
# Tests for DNA and RNA added on 20250305
info <- "read_fasta() handles 0-length 'iseq' argument for DNA"
no_DNA_dat <- read_fasta(fasta_file, iseq = numeric(), molecule = "DNA")
expect_equal(nrow(no_DNA_dat), 0, info = info)
expect_equal(colnames(no_DNA_dat), c("gene", "organism", "ref", "abbrv", "chains", "A", "C", "G", "T"), info = info)
info <- "read_fasta() handles 0-length 'iseq' argument for RNA"
no_RNA_dat <- read_fasta(fasta_file, iseq = numeric(), molecule = "RNA")
expect_equal(nrow(no_RNA_dat), 0, info = info)
expect_equal(colnames(no_RNA_dat), c("gene", "organism", "ref", "abbrv", "chains", "A", "C", "G", "U"), info = info)

# Moved from CHNOSZ on 20240328
info <- "read_fasta() reads selected sequences correctly"
fasta_file <- system.file("extdata/fasta/redoxin.fasta", package = "canprot")
aa <- read_fasta(fasta_file)
aa1 <- read_fasta(fasta_file, 1)
expect_equal(aa1, aa[1, ], info = info)
aa8rev <- read_fasta(fasta_file, 8:1)
# Use unlist here so that different row names are not compared
expect_equal(unlist(aa8rev), unlist(aa[8:1, ]), info = info)

# Tests added on 20250304

info <- "Asking for lines returns same output as readLines()"
fasta_file <- system.file("extdata/fasta/redoxin.fasta", package = "canprot")
fastalines <- read_fasta(fasta_file, type = "lines")
filelines <- readLines(fasta_file)
expect_equal(fastalines, filelines, info = info)

info <- "Asking for sequences returns sequences"
seqs <- read_fasta(fasta_file, type = "seqs")
expect_equal(range(sapply(seqs, nchar)), c(83, 508), info = info)

info <- "Reading DNA sequence produces expected message"
expect_message(DNAdat <- read_fasta(fasta_file, molecule = "DNA"), "unrecognized letter\\(s\\) in DNA sequence: D E F H I K L M N P Q R S V W Y", info = info)
info <- "Reading DNA sequence returns expected columns"
expect_equal(colnames(DNAdat)[6:9], c("A", "C", "G", "T"), info = info)
info <- "Reading RNA sequence produces expected message"
expect_message(RNAdat <- read_fasta(fasta_file, molecule = "RNA"), "unrecognized letter\\(s\\) in RNA sequence: D E F H I K L M N P Q R S T V W Y", info = info)
info <- "Reading RNA sequence returns expected columns"
expect_equal(colnames(RNAdat)[6:9], c("A", "C", "G", "U"), info = info)
info <- "Stop for unknown molecule"
expect_error(read_fasta(fasta_file, molecule = "AAA"), "unknown molecule 'AAA'", info = info)
