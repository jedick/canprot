# canprot/R/pdat_lung.R
# retrieve protein IDs for lung cancer studies
# 20160717-20190408 assemble data for 2020 compilation

pdat_lung <- function(dataset = 2020, basis = "rQEC") {
  if(identical(dataset, 2020)) {
    return(c(
             "LXC+06",
             "KNT+10",
             "KHA+12_ADC", "KHA+12_SCC", "YLL+12", "ZZD+12",
             "ZZY+13",
             "LLY+14", "LWT+14", "ZLH+14",
             "FGP+16", "JCP+16",
             "HHH+16_pN0=PP", "HHH+16_pN1=PP", "HHH+16_pN2.M1=PP",
             "TLB+16",
             "FGW+17",
             "SFS+17_LF",
             "WLC+17",
             "YCC+17_SqCC.Oncogene", "YCC+17_SqCC.TSG", "YCC+17_SqCC.Glycoprotein",
             "YCC+17_ADC.Oncogene", "YCC+17_ADC.TSG", "YCC+17_ADC.Glycoprotein"
             ))
  }
  # remove tags
  dataset <- strsplit(dataset, "=")[[1]][1]
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse="_")
  extdatadir <- system.file("extdata", package="canprot")
  datadir <- paste0(extdatadir, "/expression/lung/")
  if(study=="KNT+10") {
    # 20160717 lung stage IIIA / I, Kawamura et al., 2010
    dat <- read.csv(paste0(datadir, "KNT+10.csv.xz"), as.is=TRUE)
    description <- "lung IIIA/I"
    dat <- check_IDs(dat, "Accession.number")
    pcomp <- protcomp(dat$Accession.number, basis=basis)
    up2 <- dat$Fold.change < 1
  } else if(study=="HHH+16") {
    # 20160720 lung stages, Hu et al., 2016
    # HHH+16_pN0, HHH+16_pN1, HHH+16_pN2.M1
    dat <- read.csv(paste0(datadir, "HHH+16.csv.xz"), as.is=TRUE)
    description <- paste(stage, "/ normal")
    # get data for the selected experiment
    icol <- match(paste0(stage, ".Normal"), colnames(dat))
    # include proteins with a large expression change
    dat <- dat[dat[, icol] <= 0.66 | dat[, icol] >= 1.5, ]
    dat <- check_IDs(dat, "Uniprot.accession.No.")
    pcomp <- protcomp(dat$Uniprot.accession.No., basis=basis)
    up2 <- dat[, icol] > 1
  } else if(study=="YLL+12") {
    # 20170113 human lung squamous carcinoma, Yan et al., 2012
    dat <- read.csv(paste0(datadir, "YLL+12.csv.xz"), as.is=TRUE)
    description <- "HLSC T/N"
    dat <- check_IDs(dat, "AC")
    up2 <- dat$X115.113 > 1
    dat <- cleanup(dat, "AC", up2)
    pcomp <- protcomp(dat$AC, basis=basis)
  } else if(study=="ZLH+14") {
    # 20170113 lung adenocarcinoma, Zhang et al., 2014
    dat <- read.csv(paste0(datadir, "ZLH+14.csv.xz"), as.is=TRUE)
    description <- "lung adenocarcinoma T/N"
    dat <- check_IDs(dat, "Accession")
    pcomp <- protcomp(dat$Accession, basis=basis)
    up2 <- dat$X117.118 > 1
  } else if(study=="LXC+06") {
    # 20190317 lung squamous carcinoma, Li et al., 2006
    dat <- read.csv(paste0(datadir, "LXC+06.csv.xz"), as.is=TRUE)
    description <- "lung squamous carcinoma"
    dat <- check_IDs(dat, "Accession.no.")
    up2 <- grepl("^T", dat$Spot.no.)
    dat <- cleanup(dat, "Accession.no.", up2)
    pcomp <- protcomp(dat$Accession.no., basis=basis)
  } else if(study=="FGP+16") {
    # 20190318 lung adenocarcinoma, Fahrmann et al., 2016
    dat <- read.csv(paste0(datadir, "FGP+16.csv.xz"), as.is=TRUE)
    description <- "lung adenocarcinoma"
    up2 <- dat$Fold.Change > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis=basis)
  } else if(study=="JCP+16") {
    # 20190318 mouse tumor / normal epithelial cells, Jin et al., 2016
    dat <- read.csv(paste0(datadir, "JCP+16.csv.xz"), as.is=TRUE)
    description <- "mouse epithelial T / N"
    up2 <- dat$Ratio.TEC.NEC > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis=basis, aa_file=paste0(extdatadir, "/aa/mouse/JCP+16_aa.csv"))
  } else if(study=="FGW+17") {
    # 20190325 lung adenocarcinoma, Fahrmann et al., 2017
    dat <- read.csv(paste0(datadir, "FGW+17.csv.xz"), as.is=TRUE)
    description <- "lung adenocarcinoma"
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- dat$Fold.Change..Tumor.Control. > 1
  } else if(study=="KHA+12") {
    # 20190325 ADC or SCC versus normal, Kikuchi et al., 2012
    # KHA+12_ADC, KHA+12_SCC
    dat <- read.csv(paste0(datadir, "KHA+12.csv.xz"), as.is=TRUE)
    description <- paste(stage, "vs normal")
    # use ADC or SCC dataset
    dat <- dat[dat$Higher.in == stage | dat$Lower.in == stage, ]
    up2 <- dat$Higher.in == stage
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis=basis)
  } else if(study=="LLY+14") {
    # 20190325 SCC / normal, Lihong et al., 2014
    dat <- read.csv(paste0(datadir, "LLY+14.csv.xz"), as.is=TRUE)
    description <- "SCC / normal"
    # remove version suffixes
    dat$Protein.ID <- sapply(strsplit(dat$Protein.ID, "\\."), "[", 1)
    dat <- check_IDs(dat, "Protein.ID")
    up2 <- dat$Average.Ratio > 0
    dat <- cleanup(dat, "Protein.ID", up2)
    pcomp <- protcomp(dat$Protein.ID, basis=basis)
  } else if(study=="LWT+14") {
    # 20190325 NSCLC tumor / normal, Li et al., 2014
    dat <- read.csv(paste0(datadir, "LWT+14.csv.xz"), as.is=TRUE)
    description <- "NSCLC tumor / normal"
    dat <- check_IDs(dat, "Swissprot.ID")
    pcomp <- protcomp(dat$Swissprot.ID, basis=basis)
    up2 <- dat$logFC > 0
  } else if(study=="SFS+17") {
    # 20190326 lung tumor / normal (method comparison), Stewart et al., 2017
    # SFS+17_LF, SFS+17_TMT, SFS+17_DIA
    dat <- read.csv(paste0(datadir, "SFS+17.csv.xz"), as.is=TRUE)
    description <- paste("lung tumor / normal", stage)
    # use data for specified method
    icol <- grep(stage, colnames(dat))
    up2 <- dat[, icol] > 0
    dat <- cleanup(dat, "UniProt", up2)
    pcomp <- protcomp(dat$UniProt, basis=basis)
  } else if(study=="TLB+16") {
    # 20190402 tumor/control, Tenzer et al., 2016
    dat <- read.csv(paste0(datadir, "TLB+16.csv.xz"), as.is=TRUE)
    description <- "tumor / control"
    # ratios in table are control/cancer: up in cancer is less than 1
    up2 <- dat$fold.change < 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis=basis)
  } else if(study=="WLC+17") {
    # 20190403 NSCLC, Wang et al., 2017
    dat <- read.csv(paste0(datadir, "WLC+17.csv.xz"), as.is=TRUE)
    description <- "NSCLC"
    # use first protein ID
    dat$Entry <- sapply(strsplit(dat$Entry, ";"), "[", 1)
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- dat$Log2.ratio..N.T.. < 0
  } else if(study=="YCC+17") {
    # 20190403 SqCC and ADC, TSG/oncogenes and glycproteins, Yang et al., 2017
    # YCC+17_SqCC.Oncogene, YCC+17_SqCC.TSG, YCC+17_SqCC.Glycoprotein
    # YCC+17_ADC.Oncogene, YCC+17_ADC.TSG, YCC+17_ADC.Glycoprotein
    dat <- read.csv(paste0(datadir, "YCC+17.csv.xz"), as.is=TRUE)
    description <- stage
    # use specified dataset
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- dat[, icol] >= 1.5
  } else if(study=="ZZD+12") {
    # 20190408 SCC, Zeng et al., 2012
    dat <- read.csv(paste0(datadir, "ZZD+12.csv.xz"), as.is=TRUE)
    description <- "SCC T / N"
    dat <- check_IDs(dat, "UniProt")
    up2 <- dat$NBE.vs..LSCC > 1
    dat <- cleanup(dat, "UniProt", up2)
    pcomp <- protcomp(dat$UniProt, basis=basis)
  } else if(study=="ZZY+13") {
    # 20190408 AdC vs normal, Zhang et al., 2013
    dat <- read.csv(paste0(datadir, "ZZY+13.csv.xz"), as.is=TRUE)
    description <- "AdC vs normal"
    pcomp <- protcomp(dat$Accession.no., basis=basis)
    up2 <- dat$AdC.vs..PNLT > 1
  } else stop(paste("lung dataset", dataset, "not available"))
  print(paste0("pdat_lung: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20190407
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, pcomp=pcomp, up2=up2, description=description))
}
