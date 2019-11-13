# canprot/R/pdat_osmotic.R
# retrieve protein IDs for hyperosmotic experiments
# 20160926 jmd

pdat_osmotic <- function(dataset=NULL, basis="QEC") {
  if(is.null(dataset)) {
    return(c(
             "PW08_2h", "PW08_10h", "PW08_12h",
             "WCM+09",
             "OBBH11=ASC",
             "CCC+12_25mM", "CCC+12_100mM",
             "KKG+12_25C_aw0.985", "KKG+12_14C_aw0.985", "KKG+12_25C_aw0.967", "KKG+12_14C_aw0.967",
             "CCCC13_25mM", "CCCC13_100mM", "TSZ+13",
             "GSC14_t30a", "GSC14_t30b", "GSC14_t30c",
             "CLG+15",
             "KLB+15_trans-suc", "KLB+15_trans-NaCl", "KLB+15_prot-suc", "KLB+15_prot-NaCl",
             "LDB+15_all", "LDB+15_high", "YDZ+15",
             "RBP+16"
             ))
  }
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse="_")
  extdatadir <- system.file("extdata", package="canprot")
  datadir <- paste0(extdatadir, "/expression/osmotic/")
  if(study=="KKG+12") {
    # 20160918 Escherichia coli, Kocharunchitt et al., 2012
    # KKG+12_25C_aw0.985, KKG+12_14C_aw0.985, KKG+12_25C_aw0.967, KKG+12_14C_aw0.967
    dat <- read.csv(paste0(datadir, "KKG+12.csv.xz"), as.is=TRUE)
    description <- paste("ECO57", stage)
    print(paste0("pdat_osmotic: ", description, " [", dataset, "]"))
    # use specified temperature and subcellular fraction
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    # drop missing proteins
    dat <- remove_entries(dat, is.na(dat$Entry), dataset, "missing")
    pcomp <- protcomp(dat$Entry, basis=basis, aa_file=paste0(extdatadir, "/aa/bacteria/KKG+12_aa.csv"))
    up2 <- dat[, icol] > 0
  } else if(study=="PW08") {
    # 20160918 yeast VHG
    # PW08_2h, PW08_10h, PW08_12h
    dat <- read.csv(paste0(datadir, "PW08.csv.xz"), as.is=TRUE)
    description <- paste("S. cerevisiae VHG", stage)
    print(paste0("pdat_osmotic: ", description, " [", dataset, "]"))
    # use specified population
    if(stage=="2h") icol <- grep("Ratio.115.114", colnames(dat))
    if(stage=="10h") icol <- grep("Ratio.116.114", colnames(dat))
    if(stage=="12h") icol <- grep("Ratio.117.114", colnames(dat))
    # filter p-values and ratios
    dat <- dat[dat[, icol+1] < 0.05, ]
    dat <- dat[dat[, icol] < 0.9 | dat[, icol] > 1.1, ]
    # drop missing proteins
    dat <- remove_entries(dat, is.na(dat$Entry), dataset, "missing")
    # drop proteins with differing expression patterns by spot
    ugi <- unique(dat$Entry)
    idrop <- numeric()
    for(gi in ugi) {
      igi <- which(dat$Entry == gi)
      if(length(unique(dat[igi, icol] > 1)) > 1) idrop <- c(idrop, igi)
    }
    dat <- remove_entries(dat, 1:nrow(dat) %in% idrop, dataset, "conflicting")
    # drop missing duplicated proteins
    dat <- remove_entries(dat, duplicated(dat$Entry), dataset, "duplicated")
    pcomp <- protcomp(dat$Entry, basis=basis, aa_file=paste0(extdatadir, "/aa/fungus/PW08_aa.csv"))
    up2 <- dat[, icol] > 1
  } else if(study=="CCC+12") {
    # 20160925 ARPE-19 retinal pigmented epithelium, Chen et al., 2012
    # CCC+12_25mM, CCC+12_100mM
    dat <- read.csv(paste0(datadir, "CCC+12.csv.xz"), as.is=TRUE)
    description <- paste(paste("ARPE-19"), stage)
    print(paste0("pdat_osmotic: ", description, " [", dataset, "]"))
    # use proteins with difference in specified condition
    icol <- grep(stage, colnames(dat))
    dat <- dat[dat[, icol[2]] < 0.05, ]
    # drop proteins with differing expression patterns by spot
    ugi <- unique(dat$Swiss.prot.No.)
    idrop <- numeric()
    for(gi in ugi) {
      igi <- which(dat$Swiss.prot.No. == gi)
      if(length(unique(dat[igi, icol[1]] > 0)) > 1) idrop <- c(idrop, igi)
    }
    dat <- remove_entries(dat, 1:nrow(dat) %in% idrop, dataset, "conflicting")
    # drop duplicated proteins
    dat <- remove_entries(dat, duplicated(dat$Swiss.prot.No.), dataset, "duplicated")
    dat <- update_IDs(dat, "Swiss.prot.No.")
    pcomp <- protcomp(dat$Swiss.prot.No., basis=basis)
    up2 <- dat[, icol[1]] > 0
  } else if(study=="CCCC13") {
    # 20160925 Chang liver cells, Chen et al., 2013
    # CCCC13_25mM, CCCC13_100mM
    dat <- read.csv(paste0(datadir, "CCCC13.csv.xz"), as.is=TRUE)
    description <- paste(paste("Chang liver cells"), stage)
    print(paste0("pdat_osmotic: ", description, " [", dataset, "]"))
    # use proteins with difference in specified condition
    icol <- grep(stage, colnames(dat))
    dat <- dat[dat[, icol[2]] < 0.05, ]
    # drop unidentified proteins
    dat <- remove_entries(dat, is.na(dat$Swiss.prot.No.), dataset, "unidentified")
    # drop proteins with differing expression patterns by spot
    ugi <- unique(dat$Swiss.prot.No.)
    idrop <- numeric()
    for(gi in ugi) {
      igi <- which(dat$Swiss.prot.No. == gi)
      if(length(unique(dat[igi, icol[1]] > 0)) > 1) idrop <- c(idrop, igi)
    }
    dat <- remove_entries(dat, 1:nrow(dat) %in% idrop, dataset, "conflicting")
    # drop duplicated proteins
    dat <- remove_entries(dat, duplicated(dat$Swiss.prot.No.), dataset, "duplicated")
    dat <- update_IDs(dat, "Swiss.prot.No.")
    pcomp <- protcomp(dat$Swiss.prot.No., basis=basis)
    up2 <- dat[, icol[1]] > 0
  } else if(study=="CLG+15") {
    # 20160925 conjunctival epithelial cells, Chen et al., 2015
    dat <- read.csv(paste0(datadir, "CLG+15.csv.xz"), as.is=TRUE)
    description <- paste("IOBA-NHC")
    print(paste0("pdat_osmotic: ", description, " [", dataset, "]"))
    # use proteins that have same direction of change in both conditions
    dat <- dat[(dat$T1 > 1 & dat$T2 > 1) | (dat$T1 < 1 & dat$T2 < 1), ]
    pcomp <- protcomp(dat$accession..UniProtKB.Swiss.Prot., basis=basis)
    up2 <- dat$T1 > 1
  } else if(study=="LDB+15") {
    # 20160925 CHO cells, Liu et al., 2015
    # LDB+15_all, LDB+15_high
    dat <- read.csv(paste0(datadir, "LDB+15.csv.xz"), as.is=TRUE)
    description <- paste("CHO", stage)
    print(paste0("pdat_osmotic: ", description, " [", dataset, "]"))
    # if "high" change is specified, take only proteins with a high level of change at all time points
    if(stage == "high") dat <- dat[rowSums(dat[, 6:8] > 0.2) == 3 | rowSums(dat[, 6:8] < -0.2) == 3, ]
    # drop unidentified proteins
    dat <- remove_entries(dat, is.na(dat$Entry), dataset, "unidentified")
    pcomp <- protcomp(dat$Entry, basis=basis, aa_file=paste0(extdatadir, "/aa/mouse/LDB+15_aa.csv"))
    up2 <- dat$SOM.Cluster == "Cluster 1"
  } else if(study=="OBBH11") {
    # 20160925 adipose-derived stem cells, Oswald et al., 2011
    dat <- read.csv(paste0(datadir, "OBBH11.csv.xz"), as.is=TRUE)
    description <- "adipose-derived stem cells"
    print(paste0("pdat_osmotic: ", description, " [", dataset, "]"))
    dat <- update_IDs(dat, "Uniprot.Protein.Code")
    pcomp <- protcomp(dat$Uniprot.Protein.Code, basis=basis)
    up2 <- dat$Elucidator.Expression.Ratio..Treated.Control. > 1
  } else if(study=="YDZ+15") {
    # 20160926 Yarrowia lipolytica, Yang et al., 2015
    dat <- read.csv(paste0(datadir, "YDZ+15.csv.xz"), as.is=TRUE)
    description <- paste("Yarrowia lipolytica")
    print(paste0("pdat_osmotic: ", description, " [", dataset, "]"))
    # drop unidentified and unquantified proteins
    dat <- remove_entries(dat, is.na(dat$Accession.No.), dataset, "unidentified")
    dat <- remove_entries(dat, is.na(dat$Av..ratio..high.low.), dataset, "unquantified")
    pcomp <- protcomp(substr(dat$Accession.No., 4, 12), basis=basis, aa_file=paste0(extdatadir, "/aa/fungus/YDZ+15_aa.csv"))
    up2 <- dat$Av..ratio..high.low. > 0
  } else if(study=="WCM+09") {
    # 20160926 mouse pancreatic islets, Waanders et al., 2009
    dat <- read.csv(paste0(datadir, "WCM+09.csv.xz"), as.is=TRUE)
    description <- paste("mouse pancreatic islets")
    print(paste0("pdat_osmotic: ", description, " [", dataset, "]"))
    # use the first UniProt ID, without isoform suffix
    dat$Uniprot <- substr(dat$Uniprot, 1, 6)
    dat <- update_IDs(dat, "Uniprot")
    pcomp <- protcomp(dat$Uniprot, basis=basis, aa_file=paste0(extdatadir, "/aa/mouse/WCM+09_aa.csv"))
    up2 <- dat$X.24h.GLUCOSE.Control > 1
  } else if(study=="GSC14") {
    # 20160926 Saccharomyces cerevisiae, Giardina et al., 2014
    # GSC14_t30a, GSC14_t30b, GSC14_t30c
    dat <- read.csv(paste0(datadir, "GSC14.csv.xz"), as.is=TRUE)
    description <- paste("S. cerevisiae", stage)
    print(paste0("pdat_osmotic: ", description, " [", dataset, "]"))
    # get data for the selected experiment
    if(stage=="t30a") icol <- grep("115.", colnames(dat))
    if(stage=="t30b") icol <- grep("116.", colnames(dat))
    if(stage=="t30c") icol <- grep("114.", colnames(dat))
    # filter data to include only significantly differentially expressed proteins
    lowP <- dat[, icol[2]] < 0.05
    lowP[is.na(lowP)] <- FALSE
    dat <- dat[lowP, ]
    pcomp <- protcomp(substr(dat$Accession.., 4, 12), basis=basis, aa_file=paste0(extdatadir, "/aa/fungus/GSC14_aa.csv"))
    # proteins that have relatively higher expression ratio than the median
    up2 <- dat[, icol[1]] > median(dat[, icol[1]])
  } else if(study=="KLB+15") {
    # 20160926 Caulobacter crescentus, Kohler et al., 2015
    # KLB+15_trans-suc, KLB+15_trans-NaCl, KLB+15_prot-suc, KLB+15_prot-NaCl
    dat <- read.csv(paste0(datadir, "KLB+15.csv.xz"), as.is=TRUE)
    if(grepl("suc", stage)) osmoticum <- "succinate"
    if(grepl("NaCl", stage)) osmoticum <- "NaCl"
    if(grepl("trans", stage)) ome <- "tr."
    if(grepl("prot", stage)) ome <- "pr."
    description <- paste("CAUCR", osmoticum, ome)
    print(paste0("pdat_osmotic: ", description, " [", dataset, "]"))
    # use protein identified in given experiment
    icol <- grep(gsub("-", ".*", stage), colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    pcomp <- protcomp(dat$Entry, basis=basis, aa_file=paste0(extdatadir, "/aa/bacteria/KLB+15_aa.csv"))
    up2 <- dat[, icol] > 0
  } else if(study=="RBP+16") {
    # 20161112 Paracoccidioides lutzii, da Silva Rodrigues et al., 2016
    dat <- read.csv(paste0(datadir, "RBP+16.csv.xz"), as.is=TRUE)
    description <- "Paracoccidioides lutzii"
    print(paste0("pdat_osmotic: ", description, " [", dataset, "]"))
    # drop duplicated proteins
    dat <- remove_entries(dat, duplicated(dat$Entry), dataset, "duplicated")
    pcomp <- protcomp(dat$Entry, basis=basis, aa_file=paste0(extdatadir, "/aa/fungus/RBP+16_aa.csv"))
    up2 <- dat$Fold.change > 1
  } else if(study=="TSZ+13") {
    # 20161113 eel gill (Anguilla japonica), Tse et al., 2013
    dat <- read.csv(paste0(datadir, "TSZ+13.csv.xz"), as.is=TRUE)
    description <- "eel gill"
    print(paste0("pdat_osmotic: ", description, " [", dataset, "]"))
    # drop missing and duplicated proteins
    dat <- remove_entries(dat, is.na(dat$Entry), dataset, "missing")
    dat <- remove_entries(dat, duplicated(dat$Entry), dataset, "duplicated")
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- dat$Fold.Change..FW.SW. > 1
  } else stop(paste("osmotic dataset", dataset, "not available"))
  return(list(dataset=dataset, basis=basis, pcomp=pcomp, up2=up2, names=names, description=description))
}
