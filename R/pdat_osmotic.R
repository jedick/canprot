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
             "CCCC13_25mM=glucose", "CCCC13_100mM=glucose", "CCW+13=glucose",
             "CLG+15",
             "KLB+15_prot-suc", "KLB+15_prot-NaCl",
             "LDB+15_all=glucose", "YDZ+15=microbial",
             "RBP+16=microbial",
             "KAK+17=microbial",
             "JBG+18=microbial", "SMS+18_wt", "SMS+18_FGFR12.deficient",
             "IXA+19=glucose", "MGF+19_10=microbial", "MGF+19_20=microbial",
             "AST+20=microbial",
             "MHP+20_H9c2=glucose", "MHP+20_HEK=glucose",
             "MPR+20_3h.high.glucose=glucose", "MPR+20_12h.high.glucose=glucose", "MPR+20_3h.high.mannitol", "MPR+20_24h.high.mannitol"
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
    description <- paste("retinal pigmented epithelium in", gsub("mM", " mM", stage), "glucose vs 5.5 mM glucose (mannitol-balanced)")
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
  } else if(study=="CLG+15") {
    # 20160925 conjunctival epithelial cells, Chen et al., 2015
    dat <- read.csv(paste0(datadir, "CLG+15.csv.xz"), as.is=TRUE)
    description <- paste("human conjunctival epithelial cells in 380 or 480 mOsm vs 280 mOsm NaCl")
    # use proteins that have same direction of change in both conditions
    dat <- dat[(dat$T1 > 1 & dat$T2 > 1) | (dat$T1 < 1 & dat$T2 < 1), ]
    pcomp <- protcomp(dat$accession..UniProtKB.Swiss.Prot., basis=basis)
    up2 <- dat$T1 > 1
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
  } else if(study=="OBBH11") {
    # 20160925 adipose-derived stem cells, Oswald et al., 2011
    dat <- read.csv(paste0(datadir, "OBBH11.csv.xz"), as.is=TRUE)
    description <- "adipose-derived stem cells in 400 mOsm vs 300 mOsm NaCl"
    dat <- check_IDs(dat, "Uniprot.Protein.Code")
    pcomp <- protcomp(dat$Uniprot.Protein.Code, basis=basis)
    up2 <- dat$Elucidator.Expression.Ratio..Treated.Control. > 1
  } else if(study=="YDZ+15") {
    # 20160926 Yarrowia lipolytica, Yang et al., 2015
    dat <- read.csv(paste0(datadir, "YDZ+15.csv.xz"), as.is=TRUE)
    description <- paste("Yarrowia lipolytica in 4.21 osmol/kg vs 3.17 osmol/kg NaCl")
    up2 <- dat$Av..ratio..high.low. > 0
    dat <- cleanup(dat, "Accession.No.", up2)
    pcomp <- protcomp(substr(dat$Accession.No., 4, 12), basis=basis, aa_file=paste0(extdatadir, "/aa/yeast/YDZ+15_aa.csv.xz"))
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
  } else if(study=="KLB+15") {
    # 20160926 Caulobacter crescentus, Kohler et al., 2015
    # KLB+15_trans-suc, KLB+15_trans-NaCl, KLB+15_prot-suc, KLB+15_prot-NaCl
    dat <- read.csv(paste0(datadir, "KLB+15.csv.xz"), as.is=TRUE)
    if(grepl("suc", stage)) osmoticum <- "200 mM sucrose vs M2 minimal salts medium plus glucose"
    if(grepl("NaCl", stage)) osmoticum <- "40/50 mM NaCl vs M2 minimal salts medium plus glucose"
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
    description <- "Paracoccidioides lutzii in 0.1 M KCl vs medium with no added KCl"
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
    description <- "Lactobacillus fermentum with vs without 1.2% w/v bile salts"
    up2 <- dat$Fold.Change > 1
    dat <- cleanup(dat, "Protein.IDs", up2)
    pcomp <- protcomp(dat$Protein.IDs, basis=basis, aa_file=paste0(extdatadir, "/aa/bacteria/KAK+17_aa.csv.xz"))
  } else if(study=="MGF+19") {
    # 20200216 Staphylococcus aureus 10 and 20% NaCl, Ming et al., 2019
    # MGF+19_10, MGF+19_20
    dat <- read.csv(paste0(datadir, "MGF+19.csv.xz"), as.is=TRUE)
    description <- paste0("Staphylococcus aureus in ", stage, "% vs 0% NaCl")
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
    description <- "Candida albicans in 1 M NaCl vs medium with no added NaCl"
    up2 <- dat$log2.ratio > 0
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis = basis, aa_file = paste0(extdatadir, "/aa/yeast/JBG+18_aa.csv.xz"))
  } else if(study=="LTH+11") {
    # 20200406 Saccharomyces cerevisae 0.7 M NaCl, Lee et al., 2011
    dat <- read.csv(paste0(datadir, "LTH+11.csv.xz"), as.is=TRUE)
    description <- "S. cerevisae in 0.7 M NaCl vs control medium (YPD)"
    up2 <- dat$Regulation == "up"
    pcomp <- protcomp(dat$Entry, basis, aa_file = paste0(extdatadir, "/aa/yeast/LTH+11_aa.csv.xz"))
  } else if(study=="AST+20") {
    # 20200407 Lactobacillus fermentum, Ali et al., 2020
    dat <- read.csv(paste0(datadir, "AST+20.csv.xz"), as.is=TRUE)
    description <- "Lactobacillus fermentum with vs without 0.3%, 0.6%, 0.9%, 1.2%, 1.5% w/v bile salts"
    up2 <- dat$Folds.change > 1
    pcomp <- protcomp(dat$Protein.IDs, basis, aa_file = paste0(extdatadir, "/aa/bacteria/AST+20_aa.csv.xz"))
  } else if(study=="SMS+18") {
    # 20200407 mouse skin in low/high humidity, Seltmann et al., 2018
    # SMS+18_wt, SMS+18_FGFR12.deficient
    dat <- read.csv(paste0(datadir, "SMS+18.csv.xz"), as.is=TRUE)
    description <- paste("epidermal lysate of mice kept in low (40%) vs high (70%) humidity,", stage)
    icol <- grep(stage, colnames(dat))
    dat <- dat[abs(dat[, icol]) > 0.5, ]
    dat <- check_IDs(dat, "Accession", aa_file = paste0(extdatadir, "/aa/mouse/SMS+18_aa.csv.xz"))
    # abundance ratios are given as high humidity / low humidity, make up-expressed be low
    up2 <- dat[, icol] < 0
    dat <- cleanup(dat, "Accession", up2)
    pcomp <- protcomp(dat$Accession, basis, aa_file = paste0(extdatadir, "/aa/mouse/SMS+18_aa.csv.xz"))
  } else if(study=="MPR+20") {
    # 20200407 human aortic endothelial cells, Madonna et al., 2020
    # MPR+20_normal.glucose.insulin, MPR+20_3h.high.glucose, MPR+20_12h.high.glucose, MPR+20_3h.high.mannitol, MPR+20_24h.high.mannitol
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
    description <- "rat INS-1β cells in 27 mM vs 11 mM glucose"
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
  } else stop(paste("osmotic dataset", dataset, "not available"))
  print(paste0("pdat_osmotic: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20191120
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, basis=basis, pcomp=pcomp, up2=up2, description=description))
}
