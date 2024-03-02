# canprot/calc.metrics.R
# Calculate selected chemical metrics for proteins
# 20191027 initial version as canprot/metrics.R
# 20230704 adapted for chem16S/calc_metrics.R
# 20240302 moved to canprot
calc.metrics <- function(AAcomp, metrics = c("Zc", "nO2", "nH2O")) {

  ## Define objects used in various calculations
  # The number of C in each amino acid residue; calculated in CHNOSZ:
  # nC_AA <- sapply(makeup(info(info(aminoacids("")))$formula), "[", "C")
  # nC_AA <- nC_AA
  # names(nC_AA) <- aminoacids(3)
  nC_AA <- c(Ala = 3, Cys = 3, Asp = 4, Glu = 5, Phe = 9, Gly = 2, His = 6, 
    Ile = 6, Lys = 6, Leu = 6, Met = 5, Asn = 4, Pro = 5, Gln = 5, 
    Arg = 6, Ser = 3, Thr = 4, Val = 5, Trp = 11, Tyr = 9)
  # Identify columns with 3-letter abbreviations for the amino acids
  isAA <- tolower(colnames(AAcomp)) %in% tolower(names(nC_AA))
  iAA <- match(tolower(colnames(AAcomp)[isAA]), tolower(names(nC_AA)))

  values <- lapply(metrics, function(metric) {
  
    if(metric == "Zc") {
      Zc(AAcomp)
    } else if(metric == "nH2O") {
      nH2O(AAcomp)
    } else if(metric == "nO2") {
      nO2(AAcomp)
    } else if(metric == "GRAVY") {
      GRAVY(AAcomp)
    } else if(metric == "pI") {
      pI(AAcomp)
    } else if(metric == "MW") {
      MW(AAcomp)
    } else if(tolower(metric) %in% c("length", "plength")) {
      plength(AAcomp)
    } else if(metric %in% c("H/C", "H_C", "HC")) {
      HC(AAcomp)
    } else if(metric %in% c("N/C", "N_C", "NC")) {
      NC(AAcomp)
    } else if(metric %in% c("O/C", "O_C", "OC")) {
      OC(AAcomp)
    } else if(metric %in% c("S/C", "S_C", "SC")) {
      SC(AAcomp)
    } else stop(paste0("'", metric, "' is not an available metric"))

  })

  values <- do.call(cbind, values)
  colnames(values) <- metrics
  as.data.frame(values)

}
