# canprot/R/pdat_osmotic_bact.R
# retrieve protein IDs for hyperosmotic experiments
# 20160926 jmd
# 20200418 extract data for bacteria from pdat_osmotic.R

pdat_osmotic_bact <- function(dataset = 2020, basis = "rQEC") {
  if(identical(dataset, 2020)) {
    return(c(
             "PNWB09",
             "FTR+10",
             "HHB+12_ATCC.334", "HHB+12_DN.114.001", "HHB+12_Shirota", "HHB+12_F.19", "HHB+12_CRL.431", "HHB+12_Rosell.215",
             "KKG+12_25C_aw0.985", "KKG+12_14C_aw0.985", "KKG+12_25C_aw0.967", "KKG+12_14C_aw0.967",
             "LPK+13", "QHT+13_24.h", "QHT+13_48.h",
             "KLB+15_prot-suc", "KLB+15_prot-NaCl",
             "DSNM16_131C", "DSNM16_310F",
             "KAK+17", "LYS+17",
             "KSK+18", "LJC+18_wt", "LJC+18_mutant",
             "LWS+19", "MGF+19_10", "MGF+19_20",
             "AST+20", "GBR+20_CIRM129", "GBR+20_CIRM1025"
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
  } else if(study=="KLB+15") {
    # 20160926 Caulobacter crescentus, Kohler et al., 2015
    # KLB+15_trans-suc, KLB+15_trans-NaCl, KLB+15_prot-suc, KLB+15_prot-NaCl
    dat <- read.csv(paste0(datadir, "KLB+15.csv.xz"), as.is=TRUE)
    if(grepl("suc", stage)) osmoticum <- "200 mM sucrose vs M2 minimal salts medium"
    if(grepl("NaCl", stage)) osmoticum <- "40/50 mM NaCl vs M2 minimal salts medium"
    description <- paste("Caulobacter crescentus in", osmoticum)
    # use protein identified in given experiment
    icol <- grep(gsub("-", ".*", stage), colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    pcomp <- protcomp(dat$Entry, basis=basis, aa_file=paste0(extdatadir, "/aa/bacteria/KLB+15_aa.csv.xz"))
    up2 <- dat[, icol] > 0
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
  } else if(study=="AST+20") {
    # 20200407 Lactobacillus fermentum, Ali et al., 2020
    dat <- read.csv(paste0(datadir, "AST+20.csv.xz"), as.is=TRUE)
    description <- "Lactobacillus fermentum with vs without 0.3% to 1.5% w/v bile salts"
    up2 <- dat$Folds.change > 1
    pcomp <- protcomp(dat$Protein.IDs, basis, aa_file = paste0(extdatadir, "/aa/bacteria/AST+20_aa.csv.xz"))
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
  } else if(study=="KSK+18") {
    # 20200413 Acidihalobacter prosperus DSM 14174, Khaleque et al., 2018
    dat <- read.csv(file.path(datadir, "KSK+18.csv.xz"), as.is=TRUE)
    description <- "Acidihalobacter prosperus DSM 14174 30 g/L / 5 g/L NaCl"
    dat <- check_IDs(dat, "Protein", aa_file = file.path(extdatadir, "aa/bacteria/KSK+18_aa.csv.xz"))
    up2 <- dat$FCProteins > 1
    pcomp <- protcomp(dat$Protein, basis, aa_file = file.path(extdatadir, "aa/bacteria/KSK+18_aa.csv.xz"))
  } else if(study=="PNWB09") {
    # 20200416 Synechocystis sp. PCC6803, Pandhal et al., 2009
    dat <- read.csv(file.path(datadir, "PNWB09.csv.xz"), as.is=TRUE)
    description <- "Synechocystis sp. PCC6803 in 6% w/v NaCl vs no added salt"
    up2 <- dat$Expression == "increased"
    pcomp <- protcomp(dat$Entry, basis, aa_file = file.path(extdatadir, "aa/bacteria/PNWB09_aa.csv.xz"))
  } else if(study=="GBR+20") {
    # 20200416 Propionibacterium freudenreichii, Gaucher et al., 2020
    # GBR+20_CIRM129, GBR+20_CIRM1025
    dat <- read.csv(file.path(datadir, "GBR+20.csv.xz"), as.is=TRUE)
    description <- paste("Propionibacterium freudenreichii", stage, "in NaCl vs MMO")
    icol <- match(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    up2 <- dat[, icol] > 1
    pcomp <- protcomp(dat$Entry, basis, aa_file = file.path(extdatadir, "aa/bacteria/GBR+20_aa.csv.xz"))
  } else if(study=="HHB+12") {
    # 20200417 Lactobacillus casei in bile salt, Hamon et al., 2012
    # HHB+12_ATCC.334, HHB+12_DN.114.001, HHB+12_Shirota, HHB+12_F.19, HHB+12_CRL.431, HHB+12_Rosell.215
    dat <- read.csv(file.path(datadir, "HHB+12.csv.xz"), as.is=TRUE)
    description <- paste("Lactobacillus casei", stage, "with/without bile salt")
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    up2 <- dat[, icol] > 0
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis, aa_file = file.path(extdatadir, "aa/bacteria/HHB+12_aa.csv.xz"))
  } else stop(paste("osmotic_bact dataset", dataset, "not available"))
  print(paste0("pdat_osmotic_bact: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20191120
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, basis=basis, pcomp=pcomp, up2=up2, description=description))
}
