# canprot/R/pdat_CRC.R
# get protein data for colorectal cancer
# 20160703 jmd
# 20161011 updated with new data [LXM+16]; add =AD tag (adenoma as n2)
# 20170904 add =NT tag (normal tissue as n1)

pdat_CRC <- function(dataset=NULL, basis="QEC") {
  # list available datasets
  if(is.null(dataset)) { 
    return(c("WTK+08=NT",
             "AKP+10_CRC", "AKP+10_CIN", "AKP+10_MIN", "JKMF10", "XZC+10_I=NT", "XZC+10_II=NT", "ZYS+10=NT",
             "BPV+11_adenoma=AD=NT", "BPV+11_stage.I=NT", "BPV+11_stage.II=NT", "BPV+11_stage.III=NT", "BPV+11_stage.IV=NT",
             "JCF+11=NT", "MRK+11_AD.NC=AD=NT", "MRK+11_AC.AD", "MRK+11_AC.NC=NT",
             "KKL+12", "KYK+12=NT", "WOD+12=NT", "YLZ+12",
             "MCZ+13=NT",
             "KWA+14", "UNS+14=AD=NT", "WKP+14",
             "STK+15=NT", "WDO+15_A.N=AD=NT", "WDO+15_C.A", "WDO+15_C.N=NT",
             "LPL+16_ACP=AD=NT", "LPL+16_CIS=NT", "LPL+16_ICC=NT", "LXM+16=NT",
             "PHL+16_AD=AD=NT", "PHL+16_CIS=NT", "PHL+16_ICC=NT"))
  }
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse="_")
  extdatadir <- system.file("extdata", package="canprot")
  datadir <- paste0(extdatadir, "/expression/CRC/")
  if(study=="JKMF10") {
    # 20150520 up- and down-regulated CRC-associated proteins reported in 4 or more studies, from Jimenez et al., 2010
    description <- "serum biomarkers up / down"
    print(paste0("pdat.CRC: ", description, " [", dataset, "]"))
    dat <- read.csv(paste0(datadir, "JKMF10.csv"), as.is=TRUE)
    # get compositional features
    pcomp <- protcomp(dat$Uniprot.ID, basis=basis)
    up2 <- dat$Change=="UP"
    names <- dat$Gene.name
  } else if(study=="KWA+14") {
    # 20150908 chromatin-binding fraction, Knol et al., 2014
    dat <- read.csv(paste0(datadir, "KWA+14.csv"), as.is=TRUE)
    description <- "chromatin-binding C / A"
    print(paste0("pdat.CRC: ", description, " [", dataset, "]"))
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- dat$Fold.change > 0 | dat$Only.in == "CRC"
    up2[is.na(up2)] <- FALSE
    names <- dat$Symbol
  } else if(study=="STK+15") {
    # 20151004 CRC membrane-enriched proteome, Sethi et al., 2015
    dat <- read.csv(paste0(datadir, "STK+15.csv"), as.is=TRUE)
    description <- "membrane enriched T / N"
    print(paste0("pdat.CRC: ", description, " [", dataset, "]"))
    # drop uncharacterized proteins
    dat <- remove_entries(dat, is.na(dat$uniprot), dataset, "missing")
    # take out one that is listed as both up and down
    dat <- remove_entries(dat, dat$uniprot=="Q9NZM1", dataset, "conflicting")
    pcomp <- protcomp(dat$uniprot, basis=basis)
    up2 <- dat$invratio > 1
    names <- dat$gene
  } else if(study=="UNS+14") {
    # 20151005 epithelial cell signature, Uzozie et al., 2014
    dat <- read.csv(paste0(datadir, "UNS+14.csv"), as.is=TRUE)
    description <- "epithelial adenoma / normal"
    print(paste0("pdat.CRC: ", description, " [", dataset, "]"))
    # remove duplicated IDs
    dat <- remove_entries(dat, duplicated(dat$uniprot), dataset, "duplicated")
    pcomp <- protcomp(dat$uniprot, basis=basis)
    up2 <- dat$log2_fold > 0
    names <- dat$Gene
  } else if(study=="BPV+11") {
    # 20160414 CRC Besson et al., 2015
    # BPV+11_adenoma, BPV+11_stage.I, BPV+11_stage.II, BPV+11_stage.III, BPV+11_stage.IV
    dat <- read.csv(paste0(datadir, "BPV+11.csv"), as.is=TRUE)
    description <- paste0(gsub("\\.", " ", stage), " / normal")
    print(paste0("pdat.CRC: ", description, " [", dataset, "]"))
    # keep signifcantly changed proteins for the cancer stage
    istage <- match(tolower(stage), tolower(colnames(dat)))
    dat <- dat[dat[, istage] != 1, ]
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- dat[, istage] > 1
  } else if(study=="WDO+15") {
    # 20160414 Wisniewski et al., 2015
    # WDO+15_A.N, WDO+15_C.A, WDO+15_C.N
    dat <- read.csv(paste0(datadir, "WDO+15.csv"), as.is=TRUE)
    if(stage=="A.N") description <- "adenoma / normal"
    if(stage=="C.A") description <- "carcinoma / adenoma"
    if(stage=="C.N") description <- "carcinoma / normal"
    print(paste0("pdat.CRC: ", description, " [", dataset, "]"))
    # columns with the signficance and ratio
    isig <- grep(paste("Significant", stage, sep="."), colnames(dat), fixed=TRUE)
    irat <- grep(paste("ratio", stage, sep="."), colnames(dat), fixed=TRUE)
    # keep only the proteins marked with a significant change
    dat <- dat[dat[, isig]=="+", ]
    # remove the CON__ prefix
    dat$Majority.protein.IDs <- gsub("CON__", "", dat$Majority.protein.IDs)
    # list the known UniProt IDs and take the first (non-NA) match
    knownIDs <- check_ID(dat$Majority.protein.IDs)
    ID <- sapply(sapply(knownIDs, na.omit), "[", 1)
    dat$Majority.protein.IDs <- ID
    # drop duplicated proteins
    dat <- remove_entries(dat, duplicated(dat$Majority.protein.IDs), dataset, "duplicated")
    pcomp <- protcomp(dat$Majority.protein.IDs, basis=basis)
    up2 <- dat[, irat] > 0
  } else if(study=="WOD+12") {
    # 20160418 CRC tumor tissue, Wisniewski et al., 2012
    dat <- read.csv(paste0(datadir, "WOD+12.csv"), as.is=TRUE)
    description <- "T / N"
    print(paste0("pdat.CRC: ", description, " [", dataset, "]"))
    # drop missing UniProt IDs
    dat <- remove_entries(dat, is.na(dat$Uniprot), dataset, "missing")
    # list known UniProt IDs
    knownIDs <- check_ID(dat$Uniprot)
    # take the first (non-NA) match
    ID <- sapply(sapply(knownIDs, na.omit), "[", 1)
    dat$Uniprot <- ID
    # drop duplicated UniProt IDs
    dat <- remove_entries(dat, duplicated(dat$Uniprot), dataset, "duplicated")
    pcomp <- protcomp(dat$Uniprot, basis=basis)
    up2 <- dat$Median.Ratio.C.N > 1
  } else if(study=="JCF+11") {
    # 20160422 tumor vs normal, Jankova et al., 2011
    dat <- read.csv(paste0(datadir, "JCF+11.csv"), as.is=TRUE)
    description <- "T / N"
    print(paste0("pdat.CRC: ", description, " [", dataset, "]"))
    pcomp <- protcomp(dat$Accession, basis=basis)
    up2 <- dat$Av..Fold.Change > 0
    names <- dat$Name
  } else if(study=="XZC+10") {
    # 20160426 stage I and II vs normal, Xie et al., 2010
    # XZC+10_I, XZC+10_II
    dat <- read.csv(paste0(datadir, "XZC+10.csv"), as.is=TRUE)
    description <- paste("stage", stage, "/ normal")
    print(paste0("pdat.CRC: ", description, " [", dataset, "]"))
    # use data for the specified stage
    icol <- grep(paste0("Log2.", stage, ".N"), colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    # drop missing proteins
    dat <- remove_entries(dat, is.na(dat$UniProt), dataset, "missing")
    # list known UniProt IDs
    knownIDs <- check_ID(dat$UniProt)
    # take the first (non-NA) match
    ID <- sapply(sapply(knownIDs, na.omit), "[", 1)
    dat$UniProt <- ID
    # drop duplicated proteins
    dat <- remove_entries(dat, duplicated(dat$UniProt), dataset, "duplicated")
    pcomp <- protcomp(dat$UniProt, basis=basis)
    up2 <- dat[, icol] > 0
    names <- dat$Gene.Symbol
  } else if(study=="AKP+10") {
    # 20160427 adenoma ADE vs CRC, CIN, MIN, Albrethsen et al., 2010
    # AKP+10_CRC, AKP+10_CIN, AKP+10_MIN
    dat <- read.csv(paste0(datadir, "AKP+10.csv"), as.is=TRUE)
    description <- paste(stage, "nuclear matrix C / A")
    print(paste0("pdat.CRC: ", description, " [", dataset, "]"))
    # use the specified data set
    icol <- grep(paste0("Fold.Change.ADE.", stage), colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    # drop duplicated proteins
    dat <- remove_entries(dat, duplicated(dat$Entry), dataset, "duplicated")
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- dat[, icol] > 0
    names <- dat$Gene.Symbol
  } else if(study=="KKL+12") {
    # 20160428 poor / good prognosis, Kim et al., 2012
    dat <- read.csv(paste0(datadir, "KKL+12.csv"), as.is=TRUE)
    description <- "poor / good prognosis"
    print(paste0("pdat.CRC: ", description, " [", dataset, "]"))
    # drop missing proteins
    dat <- remove_entries(dat, is.na(dat$UniProt), dataset, "missing")
    # list known UniProt IDs
    knownIDs <- check_ID(dat$UniProt)
    # take the first (non-NA) match
    ID <- sapply(sapply(knownIDs, na.omit), "[", 1)
    dat$UniProt <- ID
    # drop duplicated proteins
    dat <- remove_entries(dat, duplicated(dat$UniProt), dataset, "duplicated")
    pcomp <- protcomp(dat$UniProt, basis=basis)
    up2 <- dat$protein.ratio..G.P. < 1
  } else if(study=="WKP+14") {
    # 20160428 tissue secretome, de Wit et al., 2014
    dat <- read.csv(paste0(datadir, "WKP+14.csv"), as.is=TRUE)
    description <- "tissue secretome T / N"
    print(paste0("pdat.CRC: ", description, " [", dataset, "]"))
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- dat$Fold.change > 0
    names <- dat$Gene.Symbol
  } else if(study=="KYK+12") {
    # 20160428 MSS-type CRC, Kang et al., 2012
    dat <- read.csv(paste0(datadir, "KYK+12.csv"), as.is=TRUE)
    description <- "MSS-type T / N"
    print(paste0("pdat.CRC: ", description, " [", dataset, "]"))
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- dat$mTRAQ.ratio..N.C.a < 0.5
    names <- dat$gene.name
  } else if(study=="ZYS+10") {
    # 20160430 microdissected T / N, Zhang et al., 2010
    dat <- read.csv(paste0(datadir, "ZYS+10.csv"), as.is=TRUE)
    description <- "microdissected T / N"
    print(paste0("pdat.CRC: ", description, " [", dataset, "]"))
    # drop missing proteins
    dat <- remove_entries(dat, is.na(dat$UniProt), dataset, "missing")
    # list known UniProt IDs
    knownIDs <- check_ID(dat$UniProt)
    # take the first (non-NA) match
    ID <- sapply(sapply(knownIDs, na.omit), "[", 1)
    dat$UniProt <- ID
    # drop duplicated proteins
    dat <- remove_entries(dat, duplicated(dat$UniProt), dataset, "duplicated")
    pcomp <- protcomp(dat$UniProt, basis=basis)
    up2 <- dat$Ratio..cancer.normal. > 1
  } else if(study=="MRK+11") {
    # 20160509 T / N, Mikula et al., 2011
    # MRK+11_AD.NC, MRK+11_AC.AD, MRK+11_AC.NC
    dat <- read.csv(paste0(datadir, "MRK+11.csv"), as.is=TRUE)
    if(stage=="AD.NC") description <- "adenoma / normal"
    if(stage=="AC.AD") description <- "adenocarcinoma / adenoma"
    if(stage=="AC.NC") description <- "adenocarcinoma / normal"
    print(paste0("pdat.CRC: ", description, " [", dataset, "]"))
    # the fold change and false discovery rate columns
    iFC <- grep(paste0(stage, ".FC"), colnames(dat))
    iFDR <- grep(paste0(stage, ".FDR"), colnames(dat))
    # apply the cutoffs
    isFC <- dat[, iFC] >= 3/2 | dat[, iFC] <= 2/3
    isFDR <- dat[, iFDR] <= 0.01
    dat <- dat[isFC & isFDR, ]
    # list known UniProt IDs
    knownIDs <- check_ID(dat$Swiss.ID)
    # take the first (non-NA) match
    ID <- sapply(sapply(knownIDs, na.omit), "[", 1)
    dat$Swiss.ID <- ID
    # drop duplicated proteins
    dat <- remove_entries(dat, duplicated(dat$Swiss.ID), dataset, "duplicated")
    pcomp <- protcomp(dat$Swiss.ID, basis=basis)
    up2 <- dat[, iFC] > 1
  } else if(study=="YLZ+12") {
    # 20160511 conditioned media T / N, Yao et al., 2012
    dat <- read.csv(paste0(datadir, "YLZ+12.csv"), as.is=TRUE)
    description <- "CM T / N"
    print(paste0("pdat.CRC: ", description, " [", dataset, "]"))
    # list known UniProt IDs
    knownIDs <- check_ID(dat$UniProt)
    # take the first (non-NA) match
    ID <- sapply(sapply(knownIDs, na.omit), "[", 1)
    dat$UniProt <- ID
    pcomp <- protcomp(dat$UniProt, basis=basis)
    up2 <- dat$Rsca > 0
  } else if(study=="WTK+08") {
    # 20160511 T / N, Watanabe et al., 2008
    dat <- read.csv(paste0(datadir, "WTK+08.csv"), as.is=TRUE)
    description <- "T / N"
    print(paste0("pdat.CRC: ", description, " [", dataset, "]"))
    # list known UniProt IDs
    knownIDs <- check_ID(dat$Accession.No.)
    # take the first (non-NA) match
    ID <- sapply(sapply(knownIDs, na.omit), "[", 1)
    dat$Accession.No. <- ID
    # drop duplicated proteins
    dat <- remove_entries(dat, duplicated(dat$Accession.No.), dataset, "duplicated")
    pcomp <- protcomp(dat$Accession.No., basis=basis)
    up2 <- dat$Average.T.N.ratio > 1
  } else if(study=="PHL+16") {
    # 20160602 AD/NC, CIS/NC, ICC/NC, Peng et al., 2016
    # PHL+16_AD, PHL+16_CIS, PHL+16_ICC
    dat <- read.csv(paste0(datadir, "PHL+16.csv"), as.is=TRUE)
    if(stage=="AD") description <- "AD / NC"
    if(stage=="CIS") description <- "CIS / NC"
    if(stage=="ICC") description <- "ICC / NC"
    print(paste0("pdat.CRC: ", description, " [", dataset, "]"))
    # AD/NC, CIS/NC, ICC/NC
    if(stage=="AD") ratio <- 2^dat$log.of.113.114 
    if(stage=="CIS") ratio <- 2^dat$log.of.115.114
    if(stage=="ICC") ratio <- 2^dat$log.of.116.114
    idiff <- ratio < 0.67 | ratio > 1.5
    dat <- dat[idiff, ]
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- ratio[idiff] > 1.5
  } else if(study=="LPL+16") {
    # 20160602 stromal ACP/NNCM, CIS/NNCM, ICC/NNCM, Li et al., 2016
    # LPL+16_ACP, LPL+16_CIS, LPL+16_ICC
    dat <- read.csv(paste0(datadir, "LPL+16.csv"), as.is=TRUE)
    if(stage=="ACP") description <- "stromal AD / NC"
    if(stage=="CIS") description <- "stromal CIS / NC"
    if(stage=="ICC") description <- "stromal ICC / NC"
    print(paste0("pdat.CRC: ", description, " [", dataset, "]"))
    icol <- grep(stage, colnames(dat))
    # keep only significantly changed proteins
    dat <- dat[dat[, icol] < 0.67 | dat[, icol] > 1.5, ]
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- dat[, icol] > 1.5
  } else if(study=="MCZ+13") {
    # 20160602 stromal T/N, Mu et al., 2013
    dat <- read.csv(paste0(datadir, "MCZ+13.csv"), as.is=TRUE)
    description <- "stromal T / N"
    print(paste0("pdat.CRC: ", description, " [", dataset, "]"))
    # drop missing proteins and duplicated proteins
    dat <- remove_entries(dat, is.na(dat$Entry), dataset, "missing")
    dat <- remove_entries(dat, duplicated(dat$Entry), dataset, "duplicated")
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- dat$CS.vs..NS > 0
  } else if(study=="LXM+16") {
    # 20160728 CRC, Liu et al., 2016
    dat <- read.csv(paste0(datadir, "LXM+16.csv"), as.is=TRUE)
    description <- "biopsy T / N"
    print(paste0("pdat_hypoxia: ", description, " [", dataset, "]"))
    # drop missing proteins
    dat <- remove_entries(dat, is.na(dat$UniProt), dataset, "missing")
    # list the known UniProt IDs and take the first (non-NA) match
    knownIDs <- check_ID(dat$UniProt)
    ID <- sapply(sapply(knownIDs, na.omit), "[", 1)
    dat$UniProt <- ID
    # drop unavailable and duplicated proteins
    dat <- remove_entries(dat, is.na(dat$UniProt), dataset, "unavailable")
    dat <- remove_entries(dat, duplicated(dat$UniProt), dataset, "duplicated")
    pcomp <- protcomp(dat$UniProt, basis=basis)
    up2 <- dat$Ratio..C.N. > 1
  } else stop(paste("CRC dataset", dataset, "not available"))
  return(list(dataset=dataset, basis=basis, description=description, pcomp=pcomp, up2=up2, names=names))
}
