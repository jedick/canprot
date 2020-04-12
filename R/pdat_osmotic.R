# canprot/R/pdat_osmotic.R
# retrieve protein IDs for hyperosmotic experiments
# 20160926 jmd

pdat_osmotic <- function(dataset = 2020, basis = "rQEC") {
  if(identical(dataset, 2020)) {
    return(c(
             "FTR+10=microbial",
             "LTH+11=microbial", "OBBH11",
             "KKG+12_25C_aw0.985=microbial", "KKG+12_14C_aw0.985=microbial", "KKG+12_25C_aw0.967=microbial", "KKG+12_14C_aw0.967=microbial",
             "LPK+13=microbial", "QHT+13_24.h=microbial", "QHT+13_48.h=microbial",
             "CLG+15", "KLB+15_prot-suc=microbial", "KLB+15_prot-NaCl=microbial", "YDZ+15=microbial",
             "DSNM16_131C=microbial", "DSNM16_310F=microbial", "RBP+16=microbial",
             "KAK+17=microbial", "LYS+17=microbial",
             "JBG+18=microbial", "LJC+18_wt=microbial", "LJC+18_mutant=microbial", "SMS+18_wt", "SMS+18_FGFR12.deficient",
             "LWS+19=microbial", "MGF+19_10=microbial", "MGF+19_20=microbial",
             "AST+20=microbial"
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
  } else if(study=="CLG+15") {
    # 20160925 conjunctival epithelial cells, Chen et al., 2015
    dat <- read.csv(paste0(datadir, "CLG+15.csv.xz"), as.is=TRUE)
    description <- paste("human conjunctival epithelial cells in 380 or 480 mOsm vs 280 mOsm NaCl")
    # use proteins that have same direction of change in both conditions
    dat <- dat[(dat$T1 > 1 & dat$T2 > 1) | (dat$T1 < 1 & dat$T2 < 1), ]
    pcomp <- protcomp(dat$accession..UniProtKB.Swiss.Prot., basis=basis)
    up2 <- dat$T1 > 1
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
  } else if(study=="KLB+15") {
    # 20160926 Caulobacter crescentus, Kohler et al., 2015
    # KLB+15_trans-suc, KLB+15_trans-NaCl, KLB+15_prot-suc, KLB+15_prot-NaCl
    dat <- read.csv(paste0(datadir, "KLB+15.csv.xz"), as.is=TRUE)
    if(grepl("suc", stage)) osmoticum <- "200 mM sucrose vs M2 minimal salts medium"
    if(grepl("NaCl", stage)) osmoticum <- "40/50 mM NaCl vs M2 minimal salts medium"
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
  } else if(study=="LJC+18") {
    # 20191102 Listeria monocytogenes membrane vesicles, Lee et al., 2018
    # LJC+18_wt, LJC+18_mutant
    dat <- read.csv(file.path(datadir, "LJC+18.csv.xz"), as.is=TRUE)
    description <- paste("Listeria monocytogenes", stage, "in 0.5 M NaCl vs control medium")
    # use selected dataset
    icol <- grep(stage, colnames(dat))
    dat <- dat[dat[, icol], ]
    pcomp <- protcomp(dat$Entry, basis=basis, aa_file=file.path(extdatadir, "aa/bacteria/LJC+18_aa.csv.xz"))
    up2 <- dat$condition=="salt"
  } else if(study=="KAK+17") {
    # 20191102 Lactobacillus fermentum NCDC 400 bile salt exposure, Kaur et al., 2017
    dat <- read.csv(paste0(datadir, "KAK+17.csv.xz"), as.is=TRUE)
    description <- "Lactobacillus fermentum with vs without 1.2% w/v bile salts"
    up2 <- dat$Fold.Change > 1
    dat <- cleanup(dat, "Protein.IDs", up2)
    pcomp <- protcomp(dat$Protein.IDs, basis=basis, aa_file=paste0(extdatadir, "/aa/bacteria/KAK+17_aa.csv.xz"))
  } else if(study=="FTR+10") {
    # 20161112 Corynebacterium glutamicum, Fränzel et al., 2010
    dat <- read.csv(file.path(datadir, "FTR+10.csv.xz"), as.is = TRUE)
    description <- "Corynebacterium glutamicum in 750 mM NaCl vs control medium"
    # exclude entries with any NA protein expression data
    dat <- dat[!rowSums(is.na(dat[, 4:6])) > 0, ]
    # use only proteins with consistent expression at 3 time points
    dat <- dat[abs(rowSums(sign(dat[, 4:6]))) == 3, ]
    pcomp <- protcomp(dat$Entry, basis = basis, aa_file = file.path(extdatadir, "aa/bacteria/FTR+10_aa.csv.xz"))
    up2 <- rowSums(sign(dat[, 4:6])) == 3
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
    description <- "Lactobacillus fermentum with vs without 0.3% to 1.5% w/v bile salts"
    up2 <- dat$Folds.change > 1
    pcomp <- protcomp(dat$Protein.IDs, basis, aa_file = paste0(extdatadir, "/aa/bacteria/AST+20_aa.csv.xz"))
  } else if(study=="SMS+18") {
    # 20200407 mouse skin in low/high humidity, Seltmann et al., 2018
    # SMS+18_wt, SMS+18_FGFR12.deficient
    dat <- read.csv(paste0(datadir, "SMS+18.csv.xz"), as.is=TRUE)
    description <- paste("epidermal lysate of mice kept in 40% vs 70% humidity,", stage)
    icol <- grep(stage, colnames(dat))
    dat <- dat[abs(dat[, icol]) > 0.5, ]
    dat <- check_IDs(dat, "Accession", aa_file = paste0(extdatadir, "/aa/mouse/SMS+18_aa.csv.xz"))
    # abundance ratios are given as high humidity / low humidity, make up-expressed be low
    up2 <- dat[, icol] < 0
    dat <- cleanup(dat, "Accession", up2)
    pcomp <- protcomp(dat$Accession, basis, aa_file = paste0(extdatadir, "/aa/mouse/SMS+18_aa.csv.xz"))
  } else if(study=="DSNM16") {
    # 20200408 Anabaena circinalis, D'Agostino et al., 2016
    # DSNM16_131C, DSNM16_310F
    dat <- read.csv(paste0(datadir, "DSNM16.csv.xz"), as.is=TRUE)
    description <- paste("Anabaena circinalis", stage, "in 0.45 mol/l NaCl vs medium without added NaCl")
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    up2 <- dat[, icol] > 0
    pcomp <- protcomp(dat$trembl, basis, aa_file = paste0(extdatadir, "/aa/bacteria/DSNM16_aa.csv.xz"))
  } else if(study=="QHT+13") {
    # 20200408 Synechocystis sp. PCC 6803, Qiao et al., 2013
    # QHT+13_24.h, QHT+13_48.h
    dat <- read.csv(paste0(datadir, "QHT+13.csv.xz"), as.is=TRUE)
    description <- paste("Synechocystis sp. PCC 6803 in 4% w/v vs 0% added NaCl for", gsub(".h", " h", stage))
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    up2 <- dat[, icol] > 1
    pcomp <- protcomp(dat$Entry, basis, aa_file = paste0(extdatadir, "/aa/bacteria/QHT+13_aa.csv.xz"))
  } else if(study %in% c("PW08", "WCM+09", "CCC+12", "CCCC13", "GSC14", "LDB+15")) {
    # 20200411 datasets from the 2017 compilation that have been moved to pdat_glucose.R
    return(pdat_glucose(dataset, basis))
  } else if(study=="LPK+13") {
    # 20200411 Lactobacillus johnsonii, Lee et al., 2013
    dat <- read.csv(paste0(datadir, "LPK+13.csv.xz"), as.is=TRUE)
    description <- "Lactobacillus johnsonii with vs without 0.1-0.3% bile salt"
    up2 <- dat$Regulation == "up"
    dat <- cleanup(dat, "UniProt", up2)
    pcomp <- protcomp(dat$UniProt, basis, aa_file = paste0(extdatadir, "/aa/bacteria/LPK+13_aa.csv.xz"))
  } else if(study=="LYS+17") {
    # 20200412 Lactobacillus salivarius LI01, Lv et al., 2017
    dat <- read.csv(paste0(datadir, "LYS+17.csv.xz"), as.is=TRUE)
    description <- "Lactobacillus salivarius LI01 with vs without 0.15% ox bile"
    up2 <- dat$log2fold_change > 0
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis, aa_file = paste0(extdatadir, "/aa/bacteria/LYS+17_aa.csv.xz"))
  } else if(study=="LWS+19") {
    # 20200412 Lactobacillus plantarum FS5-5, Li et al., 2019
    dat <- read.csv(paste0(datadir, "LWS+19.csv.xz"), as.is=TRUE)
    description <- "Lactobacillus plantarum FS5-5 in 6-8% w/v vs 0% NaCl"
    up2 <- dat$X118.113 > 1
    dat <- cleanup(dat, "No.", up2)
    pcomp <- protcomp(dat$No., basis, aa_file = paste0(extdatadir, "/aa/bacteria/LWS+19_aa.csv.xz"))
  } else stop(paste("osmotic dataset", dataset, "not available"))
  print(paste0("pdat_osmotic: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20191120
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, basis=basis, pcomp=pcomp, up2=up2, description=description))
}
