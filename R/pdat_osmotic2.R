# new datasets added 20191104
pdat_osmotic2 <- function(dataset=NULL, basis="rQEC") {
  if(is.null(dataset)) {
    return(c("LRB+09_2.6", "LRB+09_5.1",
             "FTR+10", "HMO+10_prot-membrane", "HMO+10_transcriptomics", # "HMO+10_prot-cytosol", # 14 up, 9 down
             "ZLZ+16_10", "ZLZ+16_17.5",
             "LLYL17_0", "LLYL17_3.5",
             "LJC+18_wt", "LJC+18_mutant",
             "JSP+19_LoS", "JSP+19_HiS" #, "JSP+19_LoT", "JSP+19_HiT"
             ))
  }
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse="_")
  extdatadir <- system.file("extdata", package="canprot")
  datadir <- paste0(extdatadir, "/expression/osmotic2/")
  if(study=="FTR+10") {
    # 20161112 Corynebacterium glutamicum, Fränzel et al., 2010
    dat <- read.csv(file.path(datadir, "FTR+10.csv.xz"), as.is = TRUE)
    description <- "Corynebacterium glutamicum"
    # exclude entries with any NA protein expression data
    dat <- dat[!rowSums(is.na(dat[, 4:6])) > 0, ]
    # use only proteins with consistent expression at 3 time points
    dat <- dat[abs(rowSums(sign(dat[, 4:6]))) == 3, ]
    pcomp <- protcomp(dat$Entry, basis = basis, aa_file = file.path(extdatadir, "aa/bacteria/FTR+10_aa.csv"))
    up2 <- rowSums(sign(dat[, 4:6])) == 3
  } else if(study=="ZLZ+16") {
    # 20191103 Nocardiopsis xinjiangensis, Zhang et al., 2016
    # ZLZ+16_10, ZLZ+16_17.5
    dat <- read.csv(file.path(datadir, "ZLZ+16.csv.xz"), as.is = TRUE)
    if(stage=="10") {
      description <- "Nocardiopsis xinjiangensis 10% / 6% NaCl"
      pval <- dat$p.Value..6..vs..10..
      icolratio <- grep("ratio_10_over_6", colnames(dat))
    }
    if(stage=="17.5") {
      description <- "Nocardiopsis xinjiangensis 17.5% / 10% NaCl"
      pval <- dat$p.Value..10..vs..17.5..
      icolratio <- grep("ratio_17.5_over_10", colnames(dat))
    }
    # use selected condition
    idiff <- (dat[, icolratio] > 1.3 | dat[, icolratio] < 1/1.3) & pval < 0.05
    dat <- dat[idiff, ]
    up2 <- dat[, icolratio] > 1.3
    pcomp <- protcomp(dat$Entry, basis = basis, aa_file = file.path(extdatadir, "aa/bacteria/ZLZ+16_aa.csv"))
  } else if(study=="LRB+09") {
    # 20191101 Halobacterium salinarum NaCl adjustment, Leuko et al., 2009
    # LRB+09_2.6, LRB+09_5.1
    dat <- read.csv(file.path(datadir, "LRB+09.csv.xz"), as.is=TRUE)
    if(stage=="2.6") description <- "H. salinarium 4.3 M / 2.6 M NaCl"
    if(stage=="5.1") description <- "H. salinarium 5.1 M / 4.3 M NaCl"
    # use selected dataset
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    # up-expressed proteins in high salinity (2.6 -> 4.3 M or 4.3 M -> 5.1 M)
    if(stage=="5.1") up2 <- dat[, icol] > 0
    if(stage=="2.6") up2 <- dat[, icol] < 0
    # drop missing proteins
    dat <- cleanup(dat, "Entry", dataset, up2)
    pcomp <- protcomp(dat$Entry, basis=basis, aa_file=file.path(extdatadir, "aa/archaea/LRB+09_aa.csv"))
  } else if(study=="LLYL17") {
    # 20191102 Tetragenococcus halophilus NaCl adjustment, Lin et al., 2017
    # LLYL_0, LLYL_3.5
    dat <- read.csv(file.path(datadir, "LLYL17.csv.xz"), as.is=TRUE)
    if(stage=="0") description <- "Tetragenococcus halophilus 1 M / 0 M NaCl"
    if(stage=="3.5") description <- "Tetragenococcus halophilus 3.5 M / 1 M NaCl"
    # use selected dataset
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    # up-expressed proteins in high salinity (1 M / 0 M or 3.5 M / 1 M)
    up2 <- dat[, icol] < 0
    pcomp <- protcomp(dat$UniProtKB.Entry, basis=basis, aa_file=file.path(extdatadir, "aa/bacteria/LLYL17_aa.csv"))
  } else if(study=="LJC+18") {
    # 20191102 Listeria monocytogenes membrane vesicles, Lee et al., 2018
    # LJC+18_wt, LJC+18_mutant
    dat <- read.csv(file.path(datadir, "LJC+18.csv.xz"), as.is=TRUE)
    description <- paste("Listeria monocytogenes", stage)
    # use selected dataset
    icol <- grep(stage, colnames(dat))
    dat <- dat[dat[, icol], ]
    pcomp <- protcomp(dat$Entry, basis=basis, aa_file=file.path(extdatadir, "aa/bacteria/LJC+18_aa.csv"))
    up2 <- dat$condition=="salt"
  } else if(study=="JSP+19") {
    # 20191102 Haloferax volcanii salt and temperature, Jevtić et al., 2019
    # JSP+19_LoS, JSP+19_HiS, JSP+19_LoT, JSP+19_HiT
    dat <- read.csv(file.path(datadir, "JSP+19.csv.xz"), as.is=TRUE)
    description <- paste("Haloferax volcanii", stage)
    # use selected condition
    icol <- grep(stage, colnames(dat))
    # at least two-fold, significant difference
    idiff <- abs(dat[, icol[1]]) >= 1 & dat[, icol[2]]==1
    dat <- dat[idiff, ]
    if(stage=="LoS") up2 <- dat[, icol[1]] < 0
    else up2 <- dat[, icol[1]] > 0
    # remove NA accessions
    dat <- cleanup(dat, "UniProt.Accession", dataset, up2)
    pcomp <- protcomp(dat$UniProt.Accession, basis=basis, aa_file=file.path(extdatadir, "aa/archaea/JSP+19_aa.csv"))
  } else if(study=="HMO+10") {
    # 20191102 Bacillus subtilis, Hahne et al., 2010
    # HMO+10_prot-cytosol, HMO+10_prot-membrane, HMO+10_transcriptomics
    dat <- read.csv(file.path(datadir, "HMO+10.csv.xz"), as.is=TRUE)
    description <- paste("Bacillus subtilis", stage)
    # use selected dataset
    dat <- dat[dat$Experiment==stage, ]
    # count up- or down-expression at each time point
    diff <- apply(sign(dat[, 4:7]), 1, sum)
    dat <- cbind(dat, diff = diff)
    # keep proteins that have same expression in at least 3 out of 4 time points
    dat <- dat[abs(dat$diff) >= 2, ]
    up2 <- dat$diff > 0
    # remove duplicated IDs
    dat <- cleanup(dat, "UniProtKB", dataset, up2)
    pcomp <- protcomp(dat$UniProtKB, basis=basis, aa_file=file.path(extdatadir, "aa/bacteria/HMO+10_aa.csv"))
  } else stop(paste("osmotic2 dataset", dataset, "not available"))
  print(paste0("pdat_osmotic2: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20190407
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, basis=basis, pcomp=pcomp, up2=up2, names=names, description=description))
}
