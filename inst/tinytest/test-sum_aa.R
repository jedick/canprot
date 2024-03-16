# Test added 20240303
AAcomp <-
structure(list(protein = c("O08452", "AMY", "AMYA", "BPT1", "CYC",
"LYSC"), organism = c("PYRFU", "BACSU", "PYRFU", "BOVIN", "BOVIN",
"CHICK"), ref = c("UniProt", "UniProt", "UniProt", "UniProt",
"UniProt", "UniProt"), abbrv = c("O08452", "P00691", "P49067",
"P00974", "P62894", "P00698"), chains = c(1L, 1L, 1L, 1L, 1L,
1L), Ala = c(28, 49, 26, 6, 6, 12), Cys = c(5, 1, 2, 6, 2, 8),
    Asp = c(33, 44, 35, 2, 3, 7), Glu = c(23, 23, 66, 2, 9, 2
    ), Phe = c(20, 20, 37, 4, 4, 3), Gly = c(45, 51, 44, 6, 14,
    12), His = c(12, 16, 14, 0, 3, 1), Ile = c(25, 35, 41, 2,
    6, 6), Lys = c(19, 30, 48, 4, 18, 6), Leu = c(27, 36, 59,
    2, 6, 8), Met = c(4, 10, 12, 1, 2, 2), Asn = c(21, 54, 24,
    3, 5, 14), Pro = c(20, 23, 28, 4, 4, 2), Gln = c(7, 29, 15,
    1, 3, 3), Arg = c(14, 24, 35, 6, 2, 11), Ser = c(21, 55,
    33, 1, 1, 10), Thr = c(16, 45, 12, 3, 8, 7), Val = c(31,
    32, 59, 1, 3, 6), Trp = c(26, 14, 17, 0, 1, 6), Tyr = c(37,
    28, 41, 4, 4, 3)), row.names = c(NA, 6L), class = "data.frame")

ref.sum <- structure(list(protein = "O08452", organism = "PYRFU", ref = "UniProt",
    abbrv = "O08452", chains = 6, Ala = 127, Cys = 24, Asp = 124,
    Glu = 125, Phe = 88, Gly = 172, His = 46, Ile = 115, Lys = 125,
    Leu = 138, Met = 31, Asn = 121, Pro = 81, Gln = 58, Arg = 92,
    Ser = 121, Thr = 91, Val = 132, Trp = 64, Tyr = 117), row.names = 1L, class = "data.frame")
ref.avg <- structure(list(protein = "O08452", organism = "PYRFU", ref = "UniProt", 
    abbrv = "O08452", chains = 1, Ala = 21.1666666666667, Cys = 4, 
    Asp = 20.6666666666667, Glu = 20.8333333333333, Phe = 14.6666666666667, 
    Gly = 28.6666666666667, His = 7.66666666666667, Ile = 19.1666666666667, 
    Lys = 20.8333333333333, Leu = 23, Met = 5.16666666666667, 
    Asn = 20.1666666666667, Pro = 13.5, Gln = 9.66666666666667, 
    Arg = 15.3333333333333, Ser = 20.1666666666667, Thr = 15.1666666666667, 
    Val = 22, Trp = 10.6666666666667, Tyr = 19.5), row.names = 1L, class = "data.frame")

expect_equal(sum_aa(AAcomp), ref.sum)
expect_equal(sum_aa(AAcomp, average = TRUE), ref.avg)
