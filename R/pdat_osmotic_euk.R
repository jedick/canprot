# canprot/R/pdat_osmotic_euk.R
# retrieve protein IDs for osmotic stress experiments in eukaryotes
# 20160926 jmd
# 20200418 eukaryotes extracted from pdat_osmotic.R

pdat_osmotic_euk <- function(dataset = 2020, basis = "rQEC") {
  if(identical(dataset, 2020)) {
    return(c(
             "LTH+11=yeast", "OBBH11",
             "LFY+12_C1h", "LFY+12_C8h", "LFY+12_C2p", "LFY+12_N1h", "LFY+12_N8h", "LFY+12_N2p",
             "CLG+15", "YDZ+15=yeast",
             "GAM+16_HTS", "GAM+16_HTS.Cmx", "RBP+16=yeast",
             "JBG+18=yeast", "SMS+18_wt", "SMS+18_FGFR12.deficient"
             ))
  }
  # remove tags
  dataset <- strsplit(dataset, "=")[[1]][1]
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse="_")
  extdatadir <- system.file("extdata", package="canprot")
  datadir <- paste0(extdatadir, "/expression/osmotic/")
  if(study=="CLG+15") {
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
  } else if(study=="RBP+16") {
    # 20161112 Paracoccidioides lutzii, da Silva Rodrigues et al., 2016
    dat <- read.csv(paste0(datadir, "RBP+16.csv.xz"), as.is=TRUE)
    description <- "Paracoccidioides lutzii in 0.1 M KCl vs medium with no added KCl"
    up2 <- dat$Fold.change > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis=basis, aa_file=paste0(extdatadir, "/aa/yeast/RBP+16_aa.csv.xz"))
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
  } else if(study=="LFY+12") {
    # 20200417 HEK293 cells, Li et al., 2012
    # LFY+12_C1h, LFY+12_C8h, LFY+12_C2p, LFY+12_N1h, LFY+12_N8h, LFY+12_N2p
    dat <- read.csv(file.path(datadir, "LFY+12.csv.xz"), as.is=TRUE)
    compartment <- ifelse(grepl("C", stage), "cytoplasm", "nucleus")
    time <- substr(stage, 2, 2)
    units <- ifelse(grepl("h", stage), "h", "passages")
    description <- paste(compartment, "of HEK293 cells in 500 (NaCl added) vs 300 mosmol/kg medium for", time, units)
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    dat <- check_IDs(dat, "UniProt.ID")
    up2 <- dat[, icol] > 0
    pcomp <- protcomp(dat$UniProt.ID, basis)
  } else if(study=="GAM+16") {
    # 20200418 human small airway epithelial cells, Gamboni et al., 2016
    # GAM+16_HTS, GAM+16_HTS.Cmx
    dat <- read.csv(file.path(datadir, "GAM+16.csv.xz"), as.is=TRUE)
    control <- paste0("isotonic", substr(stage, 4, 7))
    description <- paste("human small airway epithelial cells in", stage, "vs", control)
    icol <- grep(paste0(stage, ".Iso"), colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    dat <- dat[dat[, icol] > 2 | dat[, icol] < 0.5, ]
    up2 <- dat[, icol] > 2
    pcomp <- protcomp(dat$Uniprot.ID, basis)
  } else stop(paste("osmotic_euk dataset", dataset, "not available"))
  print(paste0("pdat_osmotic_euk: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20191120
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, basis=basis, pcomp=pcomp, up2=up2, description=description))
}
