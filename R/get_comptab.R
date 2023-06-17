# canprot/R/get_comptab.R
# Merge old Zc_nH2O() and CNS() functions, and add volume 20170718
# CNS: elemental abundance (C, N, S) per residue 20170124
# Zc_nH2O: plot and summarize Zc and nH2O/residue of proteins 20160706

get_comptab <- function(pdat, var1="Zc", var2="nH2O", plot.it=FALSE, mfun="median", oldstyle = FALSE, basis = getOption("basis")) {
  # Define functions for the possible variables of interest
  # Calculate metrics with canprot functions, not CHNOSZ 202010015
  Zc <- function() canprot::Zc(pdat$pcomp$aa)
  nH2O <- function() canprot::nH2O(pdat$pcomp$aa, basis = basis)
  nO2 <- function() canprot::nO2(pdat$pcomp$aa, basis = basis)
  # Calculate protein length 20201015
  #AA3 <- CHNOSZ::aminoacids(3)
  AA3 <- c("Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu",
           "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr")
  # The columns for the amino acids
  icol <- match(AA3, colnames(pdat$pcomp$aa))
  nAA <- function() rowSums(pdat$pcomp$aa[, icol])
  V0 <- function() {
    # Memoize the volumes so we don't depend on the CHNOSZ database 20200509
    #indices <- info(c(paste0("[", AA3, "]"), "[UPBB]", "[AABB]"))
    #volumes <- get("thermo", CHNOSZ::CHNOSZ)$obigt$V[indices]
    volumes <- c(26.864, 39.913, 41.116, 56.621, 88.545, 9.606, 65.753, 72.204, 
                 75.048, 74.2, 71.832, 43.826, 49.049, 60.078, 105.825, 27.042, 
                 44.03, 57.279, 110.045, 90.904, 26.296, 33.585)
    # Standard molal volumes of amino acid sidechains and protein backbone and terminal groups
    vAA <- volumes[1:20]
    vUPBB <- volumes[21]
    vAABB <- volumes[22]
    # The total volume of amino acid residues
    VAA <- rowSums(t(t(pdat$pcomp$aa[, icol]) * vAA))
    # The total volume of the backbone and terminal groups
    chains <- pdat$pcomp$aa$chains
    plength <- nAA()
    VUPBB <- (plength - chains) * vUPBB
    VAABB <- chains * vAABB
    # The per-residue volume (total volume of the protein divided by the length)
    (VAA + VUPBB + VAABB) / plength
  }
  # GRAVY and pI added 20191028
  # canprot:: is used to access the functions in the package namespace, not the ones defined here
  GRAVY <- function() canprot::GRAVY(pdat$pcomp$aa)
  pI <- function() canprot::pI(pdat$pcomp$aa)
  # MW (molecular weight) added 20200501
  MW <- function() MWAA(pdat$pcomp$aa)

  # Get the values of the variables using the functions
  val1 <- get(var1)()
  val2 <- get(var2)()
  if(plot.it) {
    # Set symbol shape and color
    # Hollow red square for up, filled blue circle for down in cancer
    col <- ifelse(pdat$up2, 2, 4)
    pch <- ifelse(pdat$up2, 0, 19)
    cex <- ifelse(pdat$up2, 1, 0.8)
    # Shuffle the order of points to mitigate overplotting effects
    i <- sample(1:length(val1))
    plot(val1[i], val2[i], xlab=cplab[[var1]], ylab=cplab[[var2]], col=col[i], pch=pch[i], cex=cex[i])
    title(pdat$description)
    mtext(pdat$dataset, side=4, cex=0.85, las=0, adj=0, line=-0.1)
  }
  # Calculate and print sample size
  val1_dn <- val1[!pdat$up2]
  val1_up <- val1[pdat$up2]
  val2_dn <- val2[!pdat$up2]
  val2_up <- val2[pdat$up2]
  message(paste0(pdat$dataset, " (", pdat$description ,"): n1 ", length(val1_dn), ", n2 ", length(val1_up)))
  # Calculate difference of means/medians, CLES, p-value
  mfun1_dn <- get(mfun)(val1_dn)
  mfun1_up <- get(mfun)(val1_up)
  mfun2_dn <- get(mfun)(val2_dn)
  mfun2_up <- get(mfun)(val2_up)
  val1.diff <- mfun1_up - mfun1_dn
  val2.diff <- mfun2_up - mfun2_dn
  out <- data.frame(dataset=pdat$dataset, description=pdat$description,
    n1=length(val1_dn), n2=length(val1_up),
    val1.median1=mfun1_dn, val1.median2=mfun1_up, val1.diff,
    val2.median1=mfun2_dn, val2.median2=mfun2_up, val2.diff, stringsAsFactors=FALSE)
  if(oldstyle) {
    val1.CLES <- 100*CLES(val1_dn, val1_up, distribution = NA)
    val2.CLES <- 100*CLES(val2_dn, val2_up, distribution = NA)
    val1.p.value <- val2.p.value <- NA
    if(!any(is.na(val1_dn)) & !any(is.na(val1_up))) val1.p.value <- stats::wilcox.test(val1_dn, val1_up)$p.value
    if(!any(is.na(val2_dn)) & !any(is.na(val2_up))) val2.p.value <- stats::wilcox.test(val2_dn, val2_up)$p.value
    # print summary messages
    nchar1 <- nchar(var1)
    nchar2 <- nchar(var2)
    start1 <- paste0(var1, substr("      MD ", nchar1, 10))
    start2 <- paste0(var2, substr("      MD ", nchar2, 10))
    message(paste0(start1, format(round(val1.diff, 3), nsmall=3, width=6),
                 ", CLES ", round(val1.CLES), "%",
                 ", p-value ", format(round(val1.p.value, 3), nsmall=3)))
    message(paste0(start2, format(round(val2.diff, 3), nsmall=3, width=6),
                 ", CLES ", round(val2.CLES), "%",
                 ", p-value ", format(round(val2.p.value, 3), nsmall=3), "\n"))
    out <- data.frame(dataset=pdat$dataset, description=pdat$description,
      n1=length(val1_dn), n2=length(val1_up),
      val1.median1=mfun1_dn, val1.median2=mfun1_up, val1.diff, val1.CLES, val1.p.value,
      val2.median1=mfun2_dn, val2.median2=mfun2_up, val2.diff, val2.CLES, val2.p.value, stringsAsFactors=FALSE)
  }
  # Convert colnames to use names of variables 
  colnames(out) <- gsub("val1", var1, colnames(out))
  colnames(out) <- gsub("val2", var2, colnames(out))
  # Convert colnames
  if(mfun == "mean") colnames(out) <- gsub("median", "mean", colnames(out))
  return(invisible(out))
}
