# canprot/R/pdat_prostate.R
# retrieve protein IDs for prostate cancer studies
# 20160420-20191212 assemble data for 2020 compilation

pdat_prostate <- function(dataset = 2020, basis = "rQEC") {
  if(identical(dataset, 2020)) {
    return(c(
             "GTR+08",
             "KPB+10_localized",
             "HZH+12_Gene=transcriptome", "HZH+12_Protein",
             "JHZ+13",
             "LCS+14",
             "IWT+16",
             "LAJ+18_PC", "LAJ+18_CRPC", "MAN+18", "TAK+18",
             "KRN+19_G1", "KRN+19_G2", "KRN+19_G3", "KRN+19_G4", "KRN+19_G5",
             "MMF+19_GS6", "MMF+19", "TOT+19", "ZYW+19_LG", "ZYW+19_HG"
             ))
  }
  # remove tags
  dataset <- strsplit(dataset, "=")[[1]][1]
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse="_")
  extdatadir <- system.file("extdata", package="canprot")
  datadir <- paste0(extdatadir, "/expression/prostate/")
  if(study=="KPB+10") {
    # 20160420 localized and metastatic prostate cancer, Khan et al., 2010
    # KPB+10_localized, KPB+10_metastatic
    dat <- read.csv(paste0(datadir, "KPB+10.csv.xz"), as.is=TRUE)
    if(stage=="localized") description <- "localized / benign"
    if(stage=="metastatic") description <- "metastatic / localized"
    #print(paste0("pdat_metastasis: ", description, " [", dataset, "]"))
    # choose localized or metastatic
    dat <- dat[dat[, stage] != "", ]
    up2 <- dat[, stage] == "up"
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis=basis)
    names <- dat$Gene.Symbol
  } else if(study=="IWT+16") {
    # 20170115 prostate cancer, Iglesias-Gato et al., 2016
    dat <- read.csv(paste0(datadir, "IWT+16.csv.xz"), as.is=TRUE)
    description <- "prostate T/N"
    dat <- check_IDs(dat, "Protein.IDs")
    pcomp <- protcomp(dat$Protein.IDs, basis=basis)
    up2 <- dat$Tmean.Cmean > 0
  } else if(study=="MAN+18") {
    # 20190319 PCa / BPH, Martiny et al., 2019
    dat <- read.csv(paste0(datadir, "MAN+18.csv.xz"), as.is=TRUE)
    description <- "PCa / BPH"
    dat <- check_IDs(dat, "ID.Protein")
    pcomp <- protcomp(dat$ID.Protein, basis=basis)
    up2 <- dat$Regulation == "Up"
  } else if(study=="LCS+14") {
    # 20190320 prostate glycoproteins, Liu et al., 2014
    dat <- read.csv(paste0(datadir, "LCS+14.csv.xz"), as.is=TRUE)
    description <- "glycoproteins"
    pcomp <- protcomp(dat$Uniprot, basis=basis)
    up2 <- dat$Cancer..Normal > 1
  } else if(study=="TAK+18") {
    # 20190321 PCa / BPH, Totten et al., 2019
    dat <- read.csv(paste0(datadir, "TAK+18.csv.xz"), as.is=TRUE)
    description <- "PCa / BPH"
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- dat$Delta > 0
  } else if(study=="TOT+19") {
    # 20190321 tissue cancer / healthy
    dat <- read.csv(paste0(datadir, "TOT+19.csv.xz"), as.is=TRUE)
    description <- "tissue cancer / healthy"
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- dat$Change == "increase"
  } else if(study=="LAJ+18") {
    # 20191202 PC / BPH and CRPC / BPH, Latonen et al., 2018
    # LAJ+18_PC, LAJ+18_CRPC
    dat <- read.csv(paste0(datadir, "LAJ+18.csv.xz"), as.is=TRUE)
    description <- paste(stage, "/ BPH")
    # get at least 2-fold differentially expressed proteins
    icol <- match(stage, colnames(dat))
    idiff <- dat[, icol] / dat$BPH > 2 | dat[, icol] / dat$BPH < 0.5
    dat <- dat[idiff, ]
    dat$Protein <- sapply(strsplit(dat$Protein, "\\|"), "[", 2)
    dat <- check_IDs(dat, "Protein")
    up2 <- dat[, icol] / dat$BPH > 2
    dat <- cleanup(dat, "Protein", up2)
    pcomp <- protcomp(dat$Protein, basis = basis)
  } else if(study=="GTR+08") {
    # 20191202 PCa / BPH, Garbis et al., 2008
    dat <- read.csv(paste0(datadir, "GTR+08.csv.xz"), as.is=TRUE)
    description <- "PCa / BPH"
    dat <- check_IDs(dat, "primary.accession.no.")
    pcomp <- protcomp(dat$primary.accession.no., basis = basis)
    up2 <- dat$mean.ratio > 1
  } else if(study=="JHZ+13") {
    # 20191202 PCa / benign, Jiang et al., 2013
    dat <- read.csv(paste0(datadir, "JHZ+13.csv.xz"), as.is=TRUE)
    description <- "PCa / benign"
    up2 <- dat$Ratio.fold. > 0
    dat <- cleanup(dat, "Protein_ID", up2)
    pcomp <- protcomp(dat$Protein_ID, basis = basis)
  } else if(study=="HZH+12") {
    # 20191202 PCa / benign (proteome and transcriptome), Han et al., 2012
    # HZH+12_Gene, HZH+12_Protein
    dat <- read.csv(paste0(datadir, "HZH+12.csv.xz"), as.is=TRUE)
    description <- paste("PCa / benign", stage)
    icol <- grep(paste0(stage, "_Fold_Change"), colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    dat <- check_IDs(dat, "Entry")
    up2 <- dat[, icol]  > 0
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis = basis)
  } else if(study=="ZYW+19") {
    # 20191202 Low-grade or high-grade vs normal, Zhou et al., 2019
    # ZYW+19_LG, ZYW+19_HG
    dat <- read.csv(paste0(datadir, "ZYW+19.csv.xz"), as.is=TRUE)
    description <- paste(stage, "/ adjacent normal")
    icol <- grep(stage, colnames(dat))
    dat <- dat[dat[, icol] != "", ]
    dat <- check_IDs(dat, "Majority.protein.IDs")
    up2 <- dat[, icol] == "Up"
    dat <- cleanup(dat, "Majority.protein.IDs", up2)
    pcomp <- protcomp(dat$Majority.protein.IDs, basis = basis)
  } else if(study=="MMF+19") {
    # 20191202 PCa / benign (GS=6 or GS=6 and GS>=8), Mantsiou et al., 2019
    # MMF+19_GS6, MMF+19
    dat <- read.csv(paste0(datadir, "MMF+19.csv.xz"), as.is=TRUE)
    description <- paste("PCa / adjacent benign", stage)
    if(stage == "GS6") icol <- grep("ratio.GS6.cancer.Vs.GS6.benign", colnames(dat)) else icol <- grep("ratio.cancer.Vs.benign", colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    up2 <- dat[, icol] > 1
    pcomp <- protcomp(dat$Accession, basis = basis)
  } else if(study=="KRN+19") {
    # 20191212 cancer grades vs benign prostatic hyperplasia, Kawahara et al., 2019
    # KRN+19_G1, KRN+19_G2, KRN+19_G3, KRN+19_G4, KRN+19_G5
    dat <- read.csv(paste0(datadir, "KRN+19.csv.xz"), as.is=TRUE)
    description <- paste(stage, "/ BPH")
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    # ratios are given as BPH / GX
    up2 <- dat[, icol] == "Down"
    pcomp <- protcomp(dat$Majority.protein.IDs, basis = basis)
  } else stop(paste("prostate dataset", dataset, "not available"))
  print(paste0("pdat_prostate: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20190429
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, pcomp=pcomp, up2=up2, names=names, description=description))
}

