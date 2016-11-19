# canprot/R/pdat_pancreatic.R
# retrieve protein IDs for pancreatic cancer studies
# 20160827 jmd

pdat_pancreatic <- function(dataset=NULL, basis="AA") {
  if(is.null(dataset)) {
    return(c("LHE+04",
             "CYD+05", "CGB+05",
             "CBP+07",
             "CTZ+09",
             "MLC+11", "PCS+11_PDAC", #"PCS+11_MCP", "PCS+11_SCP", 
             "TMW+11",
             "KBK+12",
             "KHO+13", "KPC+13_all", #"KPC+13_2-fold-signif",
             "PKB+13_AIP", "PKB+13_CP",
             "WLL+13_low", "WLL+13_high", "WLL+13a_PC_NT", "WLL+13a_PC.DM_NT.DM",
             "ZNWL13", "ISI+14",
             "KKC+16_T4=mouse", "KKC+16_T3=mouse", "KKC+16_T2=mouse", "KKC+16_T1=mouse"))
  }
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse="_")
  extdatadir <- system.file("extdata", package="canprot")
  datadir <- paste0(extdatadir, "/expression/pancreatic/")
  if(study=="KKC+16") {
    # 20160717 mouse PDAC, Kuo et al., 2016
    # KKC+16_T1, KKC+16_T2, KKC+16_T3, KKC+16_T4 (10, 5, 3.5, 2.5 weeks)
    dat <- read.csv(paste0(datadir, "KKC+16.csv"), as.is=TRUE)
    if(stage=="T1") description <- "mouse 10 w, T / N"
    if(stage=="T2") description <- "mouse 5 w, T / N"
    if(stage=="T3") description <- "mouse 3.5 w, T / N"
    if(stage=="T4") description <- "mouse 2.5 w, T / N"
    print(paste0("pdat_pancreatic: ", description, " [", dataset, "]"))
    # use only proteins differentially expressed at this stage
    icol <- grep(paste0(stage, ".N"), colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    pcomp <- protcomp(dat$Entry, basis=basis, aa_file=paste0(extdatadir, "/aa/mouse/KKC+16_aa.csv"))
    up2 <- dat[, icol] > 1
  } else if(study=="ISI+14") {
    # 20160827 PDAC, Iuga et al., 2014
    dat <- read.csv(paste0(datadir, "ISI+14.csv"), as.is=TRUE)
    description <- "T / N"
    print(paste0("pdat_pancreatic: ", description, " [", dataset, "]"))
    pcomp <- protcomp(dat$Protein.accesion.number, basis=basis)
    up2 <- dat$Simple.ratio > 1
  } else if(study=="MLC+11") {
    # 20160827 PDAC, McKinney et al., 2011
    dat <- read.csv(paste0(datadir, "MLC+11.csv"), as.is=TRUE)
    description <- "T / N"
    print(paste0("pdat_pancreatic: ", description, " [", dataset, "]"))
    # drop missing proteins
    dat <- remove_entries(dat, is.na(dat$UniProt), dataset, "missing")
    # list the known UniProt IDs and take the first (non-NA) match
    knownIDs <- check_ID(dat$UniProt)
    dat$UniProt <- sapply(sapply(knownIDs, na.omit), "[", 1)
    # drop ambiguous proteins
    up <- dat$UniProt[dat$Regulated == "up"]
    down <- dat$UniProt[dat$Regulated == "down"]
    iambi <- dat$UniProt %in% intersect(up, down)
    dat <- remove_entries(dat, iambi, dataset, "ambiguous")
    # drop duplicated proteins
    dat <- remove_entries(dat, duplicated(dat$UniProt), dataset, "duplicated")
    pcomp <- protcomp(dat$UniProt, basis=basis)
    up2 <- dat$Regulated == "up"
  } else if(study=="PCS+11") {
    # 20160828 PDAC, Pan et al., 2011
    # PCS+11_MCP, PCS+11_SCP, PCS+11_PDAC
    dat <- read.csv(paste0(datadir, "PCS+11.csv"), as.is=TRUE)
    description <- paste("FFPE,", stage, " / NL")
    print(paste0("pdat_pancreatic: ", description, " [", dataset, "]"))
    # keep proteins with reported ratio
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    # drop missing and duplicated proteins
    dat <- remove_entries(dat, is.na(dat$Entry), dataset, "missing")
    dat <- remove_entries(dat, duplicated(dat$Entry), dataset, "duplicated")
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- dat[, icol] > 1
  } else if(study=="WLL+13") {
    # 20160829 PDAC, Wang et al., 2013
    # WLL+13_low, WLL+13_high
    dat <- read.csv(paste0(datadir, "WLL+13.csv"), as.is=TRUE)
    description <- paste0(stage, "-grade T / N")
    print(paste0("pdat_pancreatic: ", description, " [", dataset, "]"))
    # which columns hold the expression data
    if(stage=="low") icol <- 4:7
    if(stage=="high") icol <- 8:11
    # 2 of the 4 ratios must be >=1.5 or <=0.667
    is.signif <- dat[, icol] >= 1.5 | dat[, icol] <= 0.667
    is.signif <- rowSums(is.signif) >= 2
    # 4 of the 4 ratios must be in the same direction
    is.up <- dat[, icol] > 1
    is.all.up <- apply(is.up, 1, all)
    is.all.down <- apply(!is.up, 1, all)
    # the final list
    is.diff <- is.signif & (is.all.down | is.all.up)
    dat <- dat[is.diff, ]
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- apply(dat[, icol] > 1, 1, all)
  } else if(study=="CGB+05") {
    # 20160829 PDAC, Crnogorac-Jurcevic et al., 2005
    dat <- read.csv(paste0(datadir, "CGB+05.csv"), as.is=TRUE)
    description <- "T / N"
    print(paste0("pdat_pancreatic: ", description, " [", dataset, "]"))
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- dat$Regulated == "up"
  } else if(study=="CTZ+09") {
    # 20160829 PDAC, Cui et al., 2009
    dat <- read.csv(paste0(datadir, "CTZ+09.csv"), as.is=TRUE)
    description <- "T / N"
    print(paste0("pdat_pancreatic: ", description, " [", dataset, "]"))
    # drop ambiguous proteins
    up <- dat$Swissprot.ID[dat$C.N > 1]
    down <- dat$Swissprot.ID[dat$C.N < 1]
    iambi <- dat$Swissprot.ID %in% intersect(up, down)
    dat <- remove_entries(dat, iambi, dataset, "ambiguous")
    # drop duplicated proteins
    dat <- remove_entries(dat, duplicated(dat$Swissprot.ID), dataset, "duplicated")
    pcomp <- protcomp(dat$Swissprot.ID, basis=basis)
    up2 <- dat$C.N > 1
  } else if(study=="KBK+12") {
    # 20160830 PDAC, Kojima et al., 2012
    dat <- read.csv(paste0(datadir, "KBK+12.csv"), as.is=TRUE)
    description <- "T / N"
    print(paste0("pdat_pancreatic: ", description, " [", dataset, "]"))
    pcomp <- protcomp(dat$Sequence.Id, basis=basis)
    up2 <- !(grepl("-", dat$Fold.Change..PDAC.Control.) | grepl("Adjacent", dat$Fold.Change..PDAC.Control.))
  } else if(study=="ZNWL13") {
    # 20160830 PDAC, Zhu et al., 2013
    dat <- read.csv(paste0(datadir, "ZNWL13.csv"), as.is=TRUE)
    description <- "T / N"
    print(paste0("pdat_pancreatic: ", description, " [", dataset, "]"))
    pcomp <- protcomp(dat$Accession, basis=basis)
    up2 <- dat$Up.Down == "Up"
  } else if(study=="KPC+13") {
    # 20160831 Kosanam et al., 2013
    # KPC+13_all, KPC+13_2-fold, KPC+13_2-fold-signif
    dat <- read.csv(paste0(datadir, "KPC+13.csv"), as.is=TRUE)
    description <- paste(stage, "T / N")
    print(paste0("pdat_pancreatic: ", description, " [", dataset, "]"))
    # drop missing proteins
    dat <- remove_entries(dat, is.na(dat$Entry), dataset, "missing")
    # use all listed proteins or at least 2-fold differentially expressed ones
    if(grepl("2-fold", stage)) dat <- dat[dat$PDAC.Benign.fold.change. >= 2 | dat$PDAC.Benign.fold.change. <= 0.5, ]
    if(grepl("signif", stage)) dat <- dat[dat$t.test.pvalue < 0.1, ]
    # drop ambiguous proteins
    up <- dat$Entry[dat$PDAC.Benign.fold.change. > 1]
    down <- dat$Entry[dat$PDAC.Benign.fold.change. < 1]
    iambi <- dat$Entry %in% intersect(up, down)
    dat <- remove_entries(dat, iambi, dataset, "ambiguous")
    # drop duplicated proteins
    dat <- remove_entries(dat, duplicated(dat$Entry), dataset, "duplicated")
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- dat$PDAC.Benign.fold.change. > 1
  } else if(study=="CYD+05") {
    # 20160907 Chen et al., 2005
    dat <- read.csv(paste0(datadir, "CYD+05.csv"), as.is=TRUE)
    description <- "T / N"
    print(paste0("pdat_pancreatic: ", description, " [", dataset, "]"))
    # drop missing proteins
    dat <- remove_entries(dat, is.na(dat$Entry), dataset, "missing")
    # drop duplicated proteins
    dat <- remove_entries(dat, duplicated(dat$Entry), dataset, "duplicated")
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- dat$Ratio..cancer.normal. > 1
  } else if(study=="CBP+07") {
    # 20160909 Chen et al., 2007
    dat <- read.csv(paste0(datadir, "CBP+07.csv"), as.is=TRUE)
    description <- "CP / N"
    print(paste0("pdat_pancreatic: ", description, " [", dataset, "]"))
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- dat$Ratio.CP.NL > 1
  } else if(study=="KHO+13") {
    # 20160910 Kawahara et al., 2013
    dat <- read.csv(paste0(datadir, "KHO+13.csv"), as.is=TRUE)
    description <- "T / N"
    print(paste0("pdat_pancreatic: ", description, " [", dataset, "]"))
    # drop missing proteins
    dat <- remove_entries(dat, is.na(dat$Entry), dataset, "missing")
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- rowMeans(dat[, 6:12]) > 1
  } else if(study=="LHE+04") {
    # 20160910 Lu et al., 2004
    dat <- read.csv(paste0(datadir, "LHE+04.csv"), as.is=TRUE)
    description <- "T / N"
    print(paste0("pdat_pancreatic: ", description, " [", dataset, "]"))
    # drop missing proteins
    dat <- remove_entries(dat, dat$Acc.no.=="", dataset, "missing")
    pcomp <- protcomp(dat$Acc.no., basis=basis)
    up2 <- dat$Higher.in == "cancer"
  } else if(study=="PKB+13") {
    # 20160910 Paulo et al., 2013
    # proteins exclusively identified in
    # autoimmune pancreatitis (AIP), chronic pancreatitis (CP), and pancreatic cancer (PC) cohorts
    # PKB+13_AIP, PKB+13_CP
    dat <- read.csv(paste0(datadir, "PKB+13.csv"), as.is=TRUE)
    description <- paste("PC /", stage)
    print(paste0("pdat_pancreatic: ", description, " [", dataset, "]"))
    # keep only proteins for the indicated comparison
    dat <- dat[dat$Cohort %in% c(stage, "PC"), ]
    pcomp <- protcomp(dat$UniProt.ID, basis=basis)
    up2 <- dat$Cohort == "PC"
  } else if(study=="TMW+11") {
    # 20160910 Turtoi et al., 2011
    dat <- read.csv(paste0(datadir, "TMW+11.csv"), as.is=TRUE)
    description <- "accessible T / N"
    print(paste0("pdat_pancreatic: ", description, " [", dataset, "]"))
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- dat$tumor - dat$normal > 1
  } else if(study=="WLL+13a") {
    # 20161110 PC +/- diabetes mellitus, Wang et al., 2013
    # WLL+13a_PC_NT, WLL+13a_PC.DM_NT.DM
    dat <- read.csv(paste0(datadir, "WLL+13a.csv"), as.is=TRUE)
    description <- stage
    print(paste0("pdat_pancreatic: ", description, " [", dataset, "]"))
    # which columns hold the expression data
    icol <- grep(stage, colnames(dat))
    # 3 of 4 experiments must have ratios >=1.5 or <=0.667
    nup <- rowSums(dat[, icol] > 3/2, na.rm=TRUE)
    ndn <- rowSums(dat[, icol] < 2/3, na.rm=TRUE)
    is.diff <- (nup >=3 & ndn==0) | (ndn >= 3 & nup==0)
    dat <- dat[is.diff, ]
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- apply(dat[, icol] > 1, 1, sum) >= 3
  } else stop(paste("pancreatic dataset", dataset, "not available"))
  return(list(dataset=dataset, basis=basis, pcomp=pcomp, up2=up2, names=names, description=description))
}
