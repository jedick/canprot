# canprot/R/pdat_glucose.R
# retrieve protein IDs for high-glucose experiments
# 20200411 glucose datasets extracted from osmotic.R

pdat_glucose <- function(dataset = 2020, basis = "rQEC") {
  if(identical(dataset, 2020)) {
    return(c(
             "PW08_2h=microbial", "PW08_10h=microbial", "PW08_12h=microbial",
             "WCM+09", "WFSL09",
             "CCC+12_25mM", "CCC+12_100mM",
             "CCCC13_25mM", "CCCC13_100mM", "CCW+13",
             "LDB+15_all",
             "IXA+19",
             "MHP+20_H9c2", "MHP+20_HEK",
             "MPR+20_3h.high.glucose", "MPR+20_24h.high.glucose", "MPR+20_3h.high.mannitol", "MPR+20_24h.high.mannitol"
## GSC14
             ))
  }
  # remove tags
  dataset <- strsplit(dataset, "=")[[1]][1]
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse="_")
  extdatadir <- system.file("extdata", package="canprot")
  datadir <- paste0(extdatadir, "/expression/glucose/")
  if(study=="PW08") {
    # 20160918 yeast VHG
    # PW08_2h, PW08_10h, PW08_12h
    dat <- read.csv(paste0(datadir, "PW08.csv.xz"), as.is=TRUE)
    description <- paste("S. cerevisiae in very high glucose (300 g/L) vs control (20 g/L) for", gsub("h", " ", stage))
    # use specified population
    if(stage=="2h") icol <- grep("Ratio.115.114", colnames(dat))
    if(stage=="10h") icol <- grep("Ratio.116.114", colnames(dat))
    if(stage=="12h") icol <- grep("Ratio.117.114", colnames(dat))
    # filter p-values and ratios
    dat <- dat[dat[, icol+1] < 0.05, ]
    dat <- dat[dat[, icol] < 0.9 | dat[, icol] > 1.1, ]
    # drop missing duplicated proteins
    up2 <- dat[, icol] > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis=basis, aa_file=paste0(extdatadir, "/aa/yeast/PW08_aa.csv.xz"))
  } else if(study=="CCC+12") {
    # 20160925 ARPE-19 retinal pigmented epithelium, Chen et al., 2012
    # CCC+12_25mM, CCC+12_100mM
    dat <- read.csv(paste0(datadir, "CCC+12.csv.xz"), as.is=TRUE)
    description <- paste("retinal pigmented epithelium in", gsub("mM", " mM", stage), "glucose vs 5.5 mM glucose")
    # use proteins with difference in specified condition
    icol <- grep(stage, colnames(dat))
    dat <- dat[dat[, icol[2]] < 0.05, ]
    dat <- check_IDs(dat, "Swiss.prot.No.")
    up2 <- dat[, icol[1]] > 0
    dat <- cleanup(dat, "Swiss.prot.No.", up2)
    pcomp <- protcomp(dat$Swiss.prot.No., basis=basis)
  } else if(study=="CCCC13") {
    # 20160925 Chang liver cells, Chen et al., 2013
    # CCCC13_25mM, CCCC13_100mM
    dat <- read.csv(paste0(datadir, "CCCC13.csv.xz"), as.is=TRUE)
    description <- paste("Chang liver cells in", gsub("mM", " mM", stage), " vs 5.5 mM glucose")
    # use proteins with difference in specified condition
    icol <- grep(stage, colnames(dat))
    dat <- dat[dat[, icol[2]] < 0.05, ]
    dat <- check_IDs(dat, "Swiss.prot.No.")
    up2 <- dat[, icol[1]] > 0
    dat <- cleanup(dat, "Swiss.prot.No.", up2)
    pcomp <- protcomp(dat$Swiss.prot.No., basis=basis)
  } else if(study=="LDB+15") {
    # 20160925 CHO cells, Liu et al., 2015
    # LDB+15_all, LDB+15_high
    dat <- read.csv(paste0(datadir, "LDB+15.csv.xz"), as.is=TRUE)
    description <- "Chinese hamster ovary cells in 15 g/L vs 5 g/L glucose"
    if(stage=="high") description <- paste(description, "(high)")
    # if "high" change is specified, take only proteins with a high level of change at all time points
    if(stage == "high") dat <- dat[rowSums(dat[, 6:8] > 0.2) == 3 | rowSums(dat[, 6:8] < -0.2) == 3, ]
    up2 <- dat$SOM.Cluster == "Cluster 1"
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis=basis, aa_file=paste0(extdatadir, "/aa/mouse/LDB+15_aa.csv.xz"))
  } else if(study=="WCM+09") {
    # 20160926 mouse pancreatic islets, Waanders et al., 2009
    dat <- read.csv(paste0(datadir, "WCM+09.csv.xz"), as.is=TRUE)
    description <- paste("mouse pancreatic islets in 16.7 mM vs 5.6 mM glucose")
    # use the first UniProt ID, without isoform suffix
    dat$Uniprot <- substr(dat$Uniprot, 1, 6)
    aa_file <- paste0(extdatadir, "/aa/mouse/WCM+09_aa.csv.xz")
    dat <- check_IDs(dat, "Uniprot", aa_file)
    pcomp <- protcomp(dat$Uniprot, basis = basis, aa_file = aa_file)
    up2 <- dat$X.24h.GLUCOSE.Control > 1
  } else if(study=="GSC14") {
    # 20160926 Saccharomyces cerevisiae, Giardina et al., 2014
    # GSC14_t30a, GSC14_t30b, GSC14_t30c
    dat <- read.csv(paste0(datadir, "GSC14.csv.xz"), as.is=TRUE)
    description <- paste("S. cerevisiae in glucose", stage)
    # get data for the selected experiment
    if(stage=="t30a") icol <- grep("115.", colnames(dat))
    if(stage=="t30b") icol <- grep("116.", colnames(dat))
    if(stage=="t30c") icol <- grep("114.", colnames(dat))
    # filter data to include only significantly differentially expressed proteins
    lowP <- dat[, icol[2]] < 0.05
    lowP[is.na(lowP)] <- FALSE
    dat <- dat[lowP, ]
    pcomp <- protcomp(substr(dat$Accession.., 4, 12), basis=basis, aa_file=paste0(extdatadir, "/aa/yeast/GSC14_aa.csv.xz"))
    # proteins that have relatively higher expression ratio than the median
    up2 <- dat[, icol[1]] > median(dat[, icol[1]])
  } else if(study=="MPR+20") {
    # 20200407 human aortic endothelial cells, Madonna et al., 2020
    # MPR+20_normal.glucose.insulin, MPR+20_3h.high.glucose, MPR+20_24h.high.glucose, MPR+20_3h.high.mannitol, MPR+20_24h.high.mannitol
    dat <- read.csv(paste0(datadir, "MPR+20.csv.xz"), as.is=TRUE)
    description <- paste("human aortic endothelial cells in", stage)
    icol <- grep(stage, colnames(dat))
    # control is normal.glucose.insulin or normal.glucose
    icontrol <- match("normal.glucose.insulin", colnames(dat))
    if(stage=="normal.glucose.insulin") icontrol <- match("normal.glucose", colnames(dat))
    # keep proteins that are present in only one of control or treatment
    dat <- dat[rowSums(dat[, c(icontrol, icol)])==1, ]
    up2 <- dat[, icol]
    pcomp <- protcomp(dat$UniProt, basis)
  } else if(study=="CCW+13") {
    # 20200407 rat INS-1β cells, Chen et al., 2013
    dat <- read.csv(paste0(datadir, "CCW+13.csv.xz"), as.is=TRUE)
    description <- "rat INS-1beta cells in 27 mM vs 11 mM glucose"
    dat <- check_IDs(dat, "Uniprot", aa_file = paste0(extdatadir, "/aa/rat/CCW+13_aa.csv.xz"))
    up2 <- dat$Average.Ratio > 1
    dat <- cleanup(dat, "Uniprot", up2)
    pcomp <- protcomp(dat$Uniprot, basis, aa_file = paste0(extdatadir, "/aa/rat/CCW+13_aa.csv.xz"))
  } else if(study=="IXA+19") {
    # 20200408 human aortal endothelial cells, Irshad et al., 2019
    dat <- read.csv(paste0(datadir, "IXA+19.csv.xz"), as.is=TRUE)
    description <- "human aortal endothelial cells in 20 mM vs 5 mM glucose"
    up2 <- dat$Fold.change > 1
    pcomp <- protcomp(dat$Entry, basis)
  } else if(study=="MHP+20") {
    # 20200408 rat H9c2 cells and human embryonic kidney cells, Meneses-Romero et al., 2020
    # MHP+20_H9c2, MHP+20_HEK
    dat <- read.csv(paste0(datadir, "MHP+20.csv.xz"), as.is=TRUE)
    if(stage=="H9c2") description <- "rat H9c2 cells in 30 mM vs 5 mM glucose"
    if(stage=="HEK") description <- "human embryonic kidney cells in 30 mM vs 5 mM glucose"
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    up2 <- dat[, icol] > 0
    pcomp <- protcomp(dat$Accession, basis, aa_file = paste0(extdatadir, "/aa/rat/MHP+20_aa.csv.xz"))
  } else if(study=="WFSL09") {
    # 20200408 bovine aortic endothelial cells, Wang et al., 2009
    dat <- read.csv(paste0(datadir, "WFSL09.csv.xz"), as.is=TRUE)
    description <- "bovine aortal endothelial cells in 22 mM vs 5 mM glucose"
    up2 <- dat$Ratio. > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis, aa_file = paste0(extdatadir, "/aa/cow/WFSL09_aa.csv.xz"))
  } else stop(paste("glucose dataset", dataset, "not available"))
  print(paste0("pdat_glucose: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20191120
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, basis=basis, pcomp=pcomp, up2=up2, description=description))
}
