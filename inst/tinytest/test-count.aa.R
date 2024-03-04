# Tests added on 2024-02-29 in canprot 1.1.2-22

info <- "Works for multiple sequences and gives message about unrecognized amino acids"
sequence <- c(seq1 = "MCMLXXVII", seq2 = "MMXXIV", seq3 = "acgt")
expect_message(AAcount <- count.aa(sequence), "unrecognized.*X", info = info)

info <- "Correct counts for amino acids"
expect_equal(AAcount$M, c(2, 2, 0), info = info)

info <- "Message about unrecognized DNA bases"
expect_message(DNAcount <- count.aa(sequence, molecule = "DNA"), "unrecognized.*I L M V X", info = info)

info <- "Correct counts for DNA (incl. lowercase letters)"
expect_equal(DNAcount$C, c(1, 0, 1), info = info)

info <- "Message about unrecognized RNA bases"
expect_message(RNAcount <- count.aa(sequence, molecule = "RNA"), "unrecognized.*I L M T V X", info = info)


# Test added on 2013-02-06 in CHNOSZ 1.0.0
info <- "count.aa() warns about unrecognized amino acids and performs substring operations"
expect_message(count.aa("ABCDEFGHIJ"), "count.aa: unrecognized letter\\(s\\) in protein sequence: B J", info = info)
myseq <- "AAAAAGGGGG"
expect_equal(count.aa(myseq, stop = 5)[, "G"], 0, info = info)
expect_equal(count.aa(myseq, start = 6)[, "A"], 0, info = info)
expect_equal(as.numeric(count.aa(myseq, start = 5, stop = 6)[, c("A", "G")]), c(1, 1), info = info)

# Test added on 2013-06-02 in CHNOSZ 1.0.3
info <- "Nucleobase sequences can be processed with count.aa()"
expect_message(dna <- count.aa("ABCDEFGHIJ", molecule = "DNA"), "count.aa: unrecognized letter\\(s\\) in DNA sequence: B D E F H I J", info = info)
expect_equal(as.numeric(dna), c(1, 1, 1, 0), info = info)

# Test added on 2018-02-17 in CHNOSZ 1.2.0
info <- "count.aa() correctly processes a longer nucleobase sequence"
seq <- "ATGTCCCGTTTCTTAGTTGCATTGGTTGCCGCACTTTTAGGAGTTGCAATTGAGATGTCCCTTCTCGTTCGCGCTCAGGGGCAGCAAACCTTGCTTTTGGCTGAAGAAAGCAAGCATTTGTCGCAATTGCGTCAACTGACTTTTGAAGGCACCAATGCCGAAGCGTATTGGTCGCCTGACGGGAAATGGTTGGTCTTTCAATCCACACGCCCACCTTACAAGGCTGACCAAATCTTCATCATGAGAGCGGATGGCTCGGGAGTTCGTGTCGTCAGCACGGGCAAAGGTCGTTGCACTTGTGCCTATTTCACGCCAGATGGCAAAGGCGTTATCTTTGCTACGACCCACCTTGCTGGACCAGAACCGCCGCAAGTGCCCAAACTGGACATTCCACGCTATGTTTGGGGCGTGTTCCCAAGTTACGAACTTTACCTGCGGCGTTTGGACACGATGGAACTTATCCGCTTGACCGATAACGAAGGCTACGACGCTGAAGCGACCATTTGCTGGAAGACTGGGCGAATTGTCTTCACAAGTTACCGCAATGGCGACCTTGACCTTTACAGCATGAAATTAGACGGCAGCGATTTGAAGCGATTGACGAAAACCATCGGCTACGAGGGCGGAGCGTTCTACTCGCCCGACGGGAAGCGGATTGTCTTCCGAGCCTATTTGCCAAAGACGCCTGACGAAATTGACGAATACAAGCGGTTGCTCCAGTTAGGCGTCATAAGCCCACCAAAGATGGAGTGGGTCGTCATGGACGCCGACGGTCGCAACATGAAGCAAATC"
counts <- data.frame(A = 190, C = 203, G = 211, T = 188)
expect_equal(as.numeric(count.aa(seq, molecule = "DNA")), as.numeric(counts), info = info)
