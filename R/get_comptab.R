# canprot/R/get_comptab.R
# merge old ZC_nH2O() and CNS() functions, and add volume 20170718
# CNS: elemental abundance (C, N, S) per residue 20170124
# ZC_nH2O: plot and summarize ZC and nH2O/residue of proteins 20160706

get_comptab <- function(pdat, var1="ZC", var2="nH2O", plot.it=FALSE, mfun="median") {
  # define functions for the possible variables of interest
  nH2O <- function() pdat$pcomp$residue.basis[, "H2O"]
  ZC <- function() pdat$pcomp$ZC
  nC <- function() pdat$pcomp$residue.formula[, "C"]
  nN <- function() pdat$pcomp$residue.formula[, "N"]
  nS <- function() pdat$pcomp$residue.formula[, "S"]
  #V0 <- function() suppressMessages(protein.obigt(pdat$pcomp$aa)$V) / pdat$pcomp$protein.length
  # longer code, but faster ...
  V0 <- function() {
    # the thermodynamic database entries for the amino acid residues
    # and the backbone and terminal groups
    AA3 <- aminoacids(3)
    indices <- info(c(paste0("[", AA3, "]"), "[UPBB]", "[AABB]"))
    volumes <- get("thermo", CHNOSZ::CHNOSZ)$obigt$V[indices]
    vAA <- volumes[1:20]
    vUPBB <- volumes[21]
    vAABB <- volumes[22]
    # the columns for the amino acids (== 6:25)
    icol <- match(AA3, colnames(pdat$pcomp$aa))
    # the total volume of amino acid residues
    VAA <- rowSums(t(t(pdat$pcomp$aa[, icol]) * vAA))
    # the protein length and total volume of the backbone and terminal groups
    pl <- rowSums(pdat$pcomp$aa[, icol])
    chains <- pdat$pcomp$aa$chains
    VUPBB <- (pl - chains) * vUPBB
    VAABB <- chains * vAABB
    # the per-residue volume (total volume of the protein divided by the length)
    (VAA + VUPBB + VAABB) / pl
  }
  nAA <- function() pdat$pcomp$protein.length
  # GRAVY and pI added 20191028
  GRAVY <- function() canprot::GRAVY(pdat$pcomp$aa)
  pI <- function() canprot::pI(pdat$pcomp$aa)
  # get the values of the variables using the functions
  val1 <- get(var1)()
  val2 <- get(var2)()
  if(plot.it) {
    # set symbol shape and color
    # hollow red square for up, filled blue circle for down in cancer
    col <- ifelse(pdat$up2, cpcol$red, cpcol$blue)
    pch <- ifelse(pdat$up2, 0, 19)
    cex <- ifelse(pdat$up2, 1, 0.8)
    # shuffle the order of points to mitigate overplotting effects
    i <- sample(1:length(val1))
    plot(val1[i], val2[i], xlab=cplab[[var1]], ylab=cplab[[var2]], col=col[i], pch=pch[i], cex=cex[i])
    title(pdat$description)
    mtext(pdat$dataset, side=4, cex=0.85, las=0, adj=0, line=-0.1)
  }
  # calculate and print sample size, difference of means/medians, CLES, p-value
  val1_dn <- val1[!pdat$up2]
  val1_up <- val1[pdat$up2]
  mfun1_dn <- get(mfun)(val1_dn)
  mfun1_up <- get(mfun)(val1_up)
  val1.diff <- mfun1_up - mfun1_dn
  val2_dn <- val2[!pdat$up2]
  val2_up <- val2[pdat$up2]
  mfun2_dn <- get(mfun)(val2_dn)
  mfun2_up <- get(mfun)(val2_up)
  val2.diff <- mfun2_up - mfun2_dn
  val1.p.value <- stats::wilcox.test(val1_dn, val1_up)$p.value
  val1.CLES <- 100*CLES(val1_dn, val1_up)
  val2.p.value <- stats::wilcox.test(val2_dn, val2_up)$p.value
  val2.CLES <- 100*CLES(val2_dn, val2_up)
  # print summary messages
  message(paste0(pdat$dataset, " (", pdat$description ,"): n1 ", length(val1_dn), ", n2 ", length(val1_up)))
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
  # convert colnames to use names of variables 
  colnames(out) <- gsub("val1", var1, colnames(out))
  colnames(out) <- gsub("val2", var2, colnames(out))
  # convert colnames
  if(mfun == "mean") colnames(out) <- gsub("median", "mean", colnames(out))
  return(invisible(out))
}
