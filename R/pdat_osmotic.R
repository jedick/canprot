# canprot/R/pdat_osmotic.R
# retrieve protein IDs for hyperosmotic experiments
# 20160926 jmd

pdat_osmotic <- function(dataset = 2020, basis = "rQEC") {
  # the only difference between 2017 and 2020 is the removal of transcriptomic datasets from KLB+15 20200326
  if(identical(dataset, 2020)) {
    return(c(
             "PW08_2h=microbial=glucose", "PW08_10h=microbial=glucose", "PW08_12h=microbial=glucose",
             "WCM+09=glucose",
             "LTH+11=microbial", "OBBH11",
             "CCC+12_25mM=glucose", "CCC+12_100mM=glucose",
             "KKG+12_25C_aw0.985=microbial", "KKG+12_14C_aw0.985=microbial", "KKG+12_25C_aw0.967=microbial", "KKG+12_14C_aw0.967=microbial",
             "CCCC13_25mM=glucose", "CCCC13_100mM=glucose",
             "GSC14_t30a=microbial=glucose", "GSC14_t30b=microbial=glucose", "GSC14_t30c=microbial=glucose",
             "CLG+15",
             "KLB+15_prot-suc", "KLB+15_prot-NaCl",
             "LDB+15_all=glucose", "YDZ+15=microbial",
             "RBP+16=microbial",
             "KAK+17=microbial",
             "JBG+18=microbial",
             "MGF+19_10=microbial", "MGF+19_20=microbial"
             ))
  }
  if(identical(dataset, 2017)) {
    return(c(
             "PW08_2h", "PW08_10h", "PW08_12h",
             "WCM+09",
             "OBBH11=ASC",
             "CCC+12_25mM", "CCC+12_100mM",
             "KKG+12_25C_aw0.985", "KKG+12_14C_aw0.985", "KKG+12_25C_aw0.967", "KKG+12_14C_aw0.967",
             "CCCC13_25mM", "CCCC13_100mM", "TSZ+13",
             "GSC14_t30a", "GSC14_t30b", "GSC14_t30c",
             "CLG+15", "KLB+15_trans-suc=transcriptome", "KLB+15_trans-NaCl=transcriptome", 
             "KLB+15_prot-suc", "KLB+15_prot-NaCl",
             "LDB+15_all", "LDB+15_high", "YDZ+15",
             "RBP+16"
             ))
  }
  # remove tags
  dataset <- strsplit(dataset, "=")[[1]][1]
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse="_")
  extdatadir <- system.file("extdata", package="canprot")
  datadir <- paste0(extdatadir, "/expression/osmotic/")
  if(study=="KKG+12") {
    # 20160918 Escherichia coli, Kocharunchitt et al., 2012
    # KKG+12_25C_aw0.985, KKG+12_14C_aw0.985, KKG+12_25C_aw0.967, KKG+12_14C_aw0.967
    dat <- read.csv(paste0(datadir, "KKG+12.csv.xz"), as.is=TRUE)
    description <- paste("E. coli in NaCl", stage)
    # use specified temperature and subcellular fraction
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    up2 <- dat[, icol] > 0
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis=basis, aa_file=paste0(extdatadir, "/aa/bacteria/KKG+12_aa.csv.xz"))
  } else if(study=="PW08") {
    # 20160918 yeast VHG
    # PW08_2h, PW08_10h, PW08_12h
    dat <- read.csv(paste0(datadir, "PW08.csv.xz"), as.is=TRUE)
    description <- paste("S. cerevisiae in glucose", stage)
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
    description <- paste("retinal pigmented epithelium in", stage, "glucose")
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
    description <- paste("Chang liver cells in", stage, "glucose")
    # use proteins with difference in specified condition
    icol <- grep(stage, colnames(dat))
    dat <- dat[dat[, icol[2]] < 0.05, ]
    dat <- check_IDs(dat, "Swiss.prot.No.")
    up2 <- dat[, icol[1]] > 0
    dat <- cleanup(dat, "Swiss.prot.No.", up2)
    pcomp <- protcomp(dat$Swiss.prot.No., basis=basis)
  } else if(study=="CLG+15") {
    # 20160925 conjunctival epithelial cells, Chen et al., 2015
    dat <- read.csv(paste0(datadir, "CLG+15.csv.xz"), as.is=TRUE)
    description <- paste("Human conjunctival epithelial cells in NaCl")
    # use proteins that have same direction of change in both conditions
    dat <- dat[(dat$T1 > 1 & dat$T2 > 1) | (dat$T1 < 1 & dat$T2 < 1), ]
    pcomp <- protcomp(dat$accession..UniProtKB.Swiss.Prot., basis=basis)
    up2 <- dat$T1 > 1
  } else if(study=="LDB+15") {
    # 20160925 CHO cells, Liu et al., 2015
    # LDB+15_all, LDB+15_high
    dat <- read.csv(paste0(datadir, "LDB+15.csv.xz"), as.is=TRUE)
    description <- "Chinese hamster ovary cells in glucose"
    if(stage=="high") description <- paste(description, "(high)")
    # if "high" change is specified, take only proteins with a high level of change at all time points
    if(stage == "high") dat <- dat[rowSums(dat[, 6:8] > 0.2) == 3 | rowSums(dat[, 6:8] < -0.2) == 3, ]
    up2 <- dat$SOM.Cluster == "Cluster 1"
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis=basis, aa_file=paste0(extdatadir, "/aa/mouse/LDB+15_aa.csv.xz"))
  } else if(study=="OBBH11") {
    # 20160925 adipose-derived stem cells, Oswald et al., 2011
    dat <- read.csv(paste0(datadir, "OBBH11.csv.xz"), as.is=TRUE)
    description <- "adipose-derived stem cells in NaCl"
    dat <- check_IDs(dat, "Uniprot.Protein.Code")
    pcomp <- protcomp(dat$Uniprot.Protein.Code, basis=basis)
    up2 <- dat$Elucidator.Expression.Ratio..Treated.Control. > 1
  } else if(study=="YDZ+15") {
    # 20160926 Yarrowia lipolytica, Yang et al., 2015
    dat <- read.csv(paste0(datadir, "YDZ+15.csv.xz"), as.is=TRUE)
    description <- paste("Yarrowia lipolytica in NaCl")
    up2 <- dat$Av..ratio..high.low. > 0
    dat <- cleanup(dat, "Accession.No.", up2)
    pcomp <- protcomp(substr(dat$Accession.No., 4, 12), basis=basis, aa_file=paste0(extdatadir, "/aa/yeast/YDZ+15_aa.csv.xz"))
  } else if(study=="WCM+09") {
    # 20160926 mouse pancreatic islets, Waanders et al., 2009
    dat <- read.csv(paste0(datadir, "WCM+09.csv.xz"), as.is=TRUE)
    description <- paste("mouse pancreatic islets in glucose")
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
  } else if(study=="KLB+15") {
    # 20160926 Caulobacter crescentus, Kohler et al., 2015
    # KLB+15_trans-suc, KLB+15_trans-NaCl, KLB+15_prot-suc, KLB+15_prot-NaCl
    dat <- read.csv(paste0(datadir, "KLB+15.csv.xz"), as.is=TRUE)
    if(grepl("suc", stage)) osmoticum <- "succinate"
    if(grepl("NaCl", stage)) osmoticum <- "NaCl"
    if(grepl("trans", stage)) ome <- "transcriptome"
    if(grepl("prot", stage)) ome <- ""
    description <- paste("Caulobacter crescentus in", osmoticum, ome)
    # use protein identified in given experiment
    icol <- grep(gsub("-", ".*", stage), colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    pcomp <- protcomp(dat$Entry, basis=basis, aa_file=paste0(extdatadir, "/aa/bacteria/KLB+15_aa.csv.xz"))
    up2 <- dat[, icol] > 0
  } else if(study=="RBP+16") {
    # 20161112 Paracoccidioides lutzii, da Silva Rodrigues et al., 2016
    dat <- read.csv(paste0(datadir, "RBP+16.csv.xz"), as.is=TRUE)
    description <- "Paracoccidioides lutzii in KCl"
    up2 <- dat$Fold.change > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis=basis, aa_file=paste0(extdatadir, "/aa/yeast/RBP+16_aa.csv.xz"))
  } else if(study=="TSZ+13") {
    # 20161113 eel gill (Anguilla japonica), Tse et al., 2013
    dat <- read.csv(paste0(datadir, "TSZ+13.csv.xz"), as.is=TRUE)
    description <- "eel gill"
    up2 <- dat$Fold.Change..FW.SW. > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis=basis)
  } else if(study=="KAK+17") {
    # 20191102 Lactobacillus fermentum NCDC 400 bile salt exposure, Kaur et al., 2017
    dat <- read.csv(paste0(datadir, "KAK+17.csv.xz"), as.is=TRUE)
    description <- "Lactobacillus fermentum in bile salts"
    up2 <- dat$Fold.Change > 1
    dat <- cleanup(dat, "Protein.IDs", up2)
    pcomp <- protcomp(dat$Protein.IDs, basis=basis, aa_file=paste0(extdatadir, "/aa/bacteria/KAK+17_aa.csv.xz"))
  } else if(study=="MGF+19") {
    # 20200216 Staphylococcus aureus 10 and 20% NaCl, Ming et al., 2019
    # MGF+19_10, MGF+19_20
    dat <- read.csv(paste0(datadir, "MGF+19.csv.xz"), as.is=TRUE)
    description <- paste0("Staphylococcus aureus in ", stage, "% NaCl")
    icol <- grep(stage, colnames(dat))
    # keep proteins with differential expression in the selected experiment
    idiff <- dat[, icol] > 2 | dat[, icol] < 0.5
    idiff[is.na(idiff)] <- FALSE
    dat <- dat[idiff, ]
    up2 <- dat[, icol] > 2
    pcomp <- protcomp(dat$Entry, basis=basis, aa_file = paste0(extdatadir, "/aa/bacteria/MGF+19_aa.csv.xz"))
  } else if(study=="JBG+18") {
    # 20200406 Candida albicans 1 M NaCl, Jacobsen et al., 2018
    dat <- read.csv(paste0(datadir, "JBG+18.csv.xz"), as.is=TRUE)
    description <- "Candida albicans in 1 M NaCl"
    up2 <- dat$log2.ratio > 0
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis = basis, aa_file = paste0(extdatadir, "/aa/yeast/JBG+18_aa.csv.xz"))
  } else if(study=="LTH+11") {
    # 20200406 Saccharomyces cerevisae 0.7 M NaCl, Lee et al., 2011
    dat <- read.csv(paste0(datadir, "LTH+11.csv.xz"), as.is=TRUE)
    description <- "S. cerevisae in 0.7 M NaCl"
    up2 <- dat$Regulation == "up"
    pcomp <- protcomp(dat$Entry, basis, aa_file = paste0(extdatadir, "/aa/yeast/LTH+11_aa.csv.xz"))
  } else stop(paste("osmotic dataset", dataset, "not available"))
  print(paste0("pdat_osmotic: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20191120
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, basis=basis, pcomp=pcomp, up2=up2, description=description))
}
