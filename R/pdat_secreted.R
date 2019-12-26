# canprot/R/pdat_secreted.R
# retrieve IDs for proteins secreted in hypoxia
# 20190325 extracted from pdat_hypoxia.R

pdat_secreted <- function(dataset = 2020, basis = "rQEC") {
  if(identical(dataset, 2020)) {
    return(c("BRA+10", "PTD+10_Hx48", "PTD+10_Hx72", #"PTD+10_ReOx=ReOx",
             "JVC+12",
             "KCW+13=transcriptome", "SKA+13", "SRS+13a_3", "SRS+13a_8",
             "LRS+14_Hy", # "LRS+14_ReoX",
             "YKK+14_soluble", "YKK+14_exosome",
             "RSE+16",
             "CGH+17_exosomes", "CGH+17_secretome",
             "CLY+18_secretome", "DWW+18", "FPR+18",
             "KAN+19_secretome", "NJVS19_CAM", "NJVS19_NTM", "PDT+19"))
  }
  # remove tags
  dataset <- strsplit(dataset, "=")[[1]][1]
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse="_")
  extdatadir <- system.file("extdata", package="canprot")
  datadir <- paste0(extdatadir, "/expression/secreted/")
  if(study=="LRS+14") {
    # 20160717 rat heart myoblast secretome, Li et al., 2014
    # LRS+14_Hy, LRS+14_Re
    dat <- read.csv(paste0(datadir, "LRS+14.csv.xz"), as.is=TRUE)
    if(stage=="Hy") description <- "myoblast secretome"
    if(stage=="Re") description <- "myoblast secretome reoxygenation / normoxia"
    # select proteins with differential expression in Hy or Re
    icol <- grep(paste0(stage, ".Ctrl_iTRAQ"), colnames(dat))
    up2 <- dat[, icol] > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis=basis, aa_file=paste0(extdatadir, "/aa/rat/LRS+14_aa.csv"))
  } else if(study=="RSE+16") {
    # 20160729 adipose-derived stem cells, Riis et al., 2016
    dat <- read.csv(paste0(datadir, "RSE+16.csv.xz"), as.is=TRUE)
    description <- "adipose-derived SC"
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- dat$Regulated == "up"
  } else if(study=="PTD+10") {
    # 20160801 A431 hypoxic / reoxygenated, Park et al., 2010
    # PTD+10_Hx48, PTD+10_Hx72, PTD+10_ReOx
    dat <- read.csv(paste0(datadir, "PTD+10.csv.xz"), as.is=TRUE)
    description <- paste("A431 cells", stage)
    if(stage=="Hx48") icol <- grep("115", colnames(dat))
    if(stage=="Hx72") icol <- grep("117", colnames(dat))
    if(stage=="ReOx") icol <- grep("116", colnames(dat))
    # filter ratio, p-value, EF value
    irat <- dat[, icol[1]] > 1.3 | dat[, icol[1]] < 1/1.3
    ipval <- dat[, icol[2]] < 0.05
    ief <- dat[, icol[3]] < 2
    dat <- dat[irat & ipval & ief, ]
    # drop missing entries
    up2 <- dat[, icol[1]] > 1.3
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis=basis)
  } else if(study=="BRA+10") {
    # 20160805 placental tissue secretome, Blankley et al., 2010
    dat <- read.csv(paste0(datadir, "BRA+10.csv.xz"), as.is=TRUE)
    description <- "placental secretome"
    dat$UniProt.accession <- sapply(strsplit(dat$UniProt.accession, "|", fixed=TRUE), "[", 2)
    dat <- check_IDs(dat, "UniProt.accession")
    pcomp <- protcomp(dat$UniProt.accession, basis=basis)
    up2 <- dat$Fold.change > 0
  } else if(study=="DWW+18") {
    # 20190322 hypoxia-induced exosomes, Dorayappan et al., 2018
    dat <- read.csv(paste0(datadir, "DWW+18.csv.xz"), as.is=TRUE)
    description <- "hypoxia-induced exosomes"
    dat <- check_IDs(dat, "Genes.symbol")
    pcomp <- protcomp(dat$Genes.symbol, basis=basis)
    up2 <- dat$FC > 1
  } else if(study=="CGH+17") {
    # 20190324 mouse cardiac fibroblast exosomes, secretome, Cosme et al., 2017
    # CGH+17_exosomes, CGH+17_secretome
    pdat <- pdat_multi(dataset, basis)
    pcomp <- pdat$pcomp
    up2 <- pdat$up2
    description <- pdat$description
    dat <- NULL
  } else if(study=="CLY+18") {
    # 20190324 HCT116 cells, Chen et al., 2018
    # CLY+18_secretome
    pdat <- pdat_multi(dataset, basis)
    pcomp <- pdat$pcomp
    up2 <- pdat$up2
    description <- pdat$description
    dat <- NULL
  } else if(study=="PDT+19") {
    # 20190326 tumor exosomes, Park et al., 2019
    dat <- read.csv(paste0(datadir, "PDT+19.csv.xz"), as.is=TRUE)
    description <- "mouse melanoma B16-F0 exosomes"
    # drop isoform suffixes
    dat$Accession <- sapply(strsplit(dat$Accession, "-"), "[", 1)
    dat <- check_IDs(dat, "Accession", aa_file=paste0(extdatadir, "/aa/mouse/PDT+19_aa.csv"))
    up2 <- dat$Log2.127.126. > 0
    dat <- cleanup(dat, "Accession", up2)
    pcomp <- protcomp(dat$Accession, basis=basis, aa_file=paste0(extdatadir, "/aa/mouse/PDT+19_aa.csv"))
  } else if(study=="SRS+13a") {
    # 20190327 placental mesenchymal stem cells 3% and 8% vs 1% O2, Salomon et al., 2013
    # SRS+13a_3, SRS+13a_8
    dat <- read.csv(paste0(datadir, "SRS+13a.csv.xz"), as.is=TRUE)
    description <- paste("pMSC", stage, "/ 1 % O2")
    # use selected dataset
    dat <- dat[dat$O2.percent %in% c(1, stage), ]
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- dat$O2.percent == stage
  } else if(study=="JVC+12") {
    # 20191204 endothelial cell-derived exosomes, de Jong et al., 2012
    dat <- read.csv(paste0(datadir, "JVC+12.csv.xz"), as.is=TRUE)
    description <- "endothelial cell-derived exosomes"
    # keep highly differential proteins
    dat <- dat[abs(dat$Hypoxia.median) > 0.2, ]
    up2 <- dat$Hypoxia.median > 0
    pcomp <- protcomp(dat$Entry, basis=basis)
  } else if(study=="KCW+13") {
    # 20191206 glioma cells, gene expression, Kucharzewska et al., 2013
    dat <- read.csv(paste0(datadir, "KCW+13.csv.xz"), as.is=TRUE)
    description <- "glioma cells gene expression"
    up2 <- dat$Fold.change > 1
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis=basis)
  } else if(study=="SKA+13") {
    # 20191207 cytotrophoblast-derived exosomes, Salomon et al., 2013
    dat <- read.csv(paste0(datadir, "SKA+13.csv.xz"), as.is=TRUE)
    description <- "cytotrophoblast-derived exosomes"
    # compare 8% / 1% conditions
    dat <- dat[dat$O2.1_percent | dat$O2.8_percent, ]
    dat <- dat[xor(dat$O2.1_percent, dat$O2.8_percent), ]
    up2 <- dat$O2.8_percent
    pcomp <- protcomp(dat$Entry, basis=basis)
  } else if(study=="YKK+14") {
    # 20191207 U373MG glioma cells, Yoon et al., 2014
    # YKK+14_soluble, YKK+14_exosome
    dat <- read.csv(paste0(datadir, "YKK+14.csv.xz"), as.is=TRUE)
    description <- paste("U373MG cells", stage)
    # get differential proteins for specified condition
    icol <- grep(stage, colnames(dat))
    dat <- dat[abs(dat[, icol]) > 0.5, ]
    up2 <- dat[, icol] > 0.5
    dat <- cleanup(dat, "Uniprot.Acc", up2)
    pcomp <- protcomp(dat$Uniprot.Acc, basis=basis)
  } else if(study=="NJVS19") {
    # 20191226 cancer-associated and normal tissue myofibroblasts, Najgebauer et al., 2019
    # NJVS19_CAM, NJVS19_NTM
    dat <- read.csv(paste0(datadir, "NJVS19.csv.xz"), as.is=TRUE)
    if(stage=="CAM") description <- "cancer-associated myofibroblasts"
    if(stage=="NTM") description <- "normal tissue myofibroblasts"
    # use selected dataset
    dat <- dat[!is.na(dat[, stage]), ]
    dat <- check_IDs(dat, "Majority.protein.IDs")
    pcomp <- protcomp(dat$Majority.protein.IDs, basis=basis)
    up2 <- dat[, stage] > 1
  } else if(study=="FPR+18") {
    # 20191226 endothelial progenitor cells, Felice et al., 2018
    dat <- read.csv(paste0(datadir, "FPR+18.csv.xz"), as.is=TRUE)
    description <- "endothelial progenitor cells"
    up2 <- dat$Modulation == "UP"
    pcomp <- protcomp(dat$Accession.number, basis=basis)
  } else if(study=="KAN+19") {
    # 20191226 human umbilical vein ECs, Kugeratski et al., 2019
    # KAN+19_secretome
    pdat <- pdat_multi(dataset, basis)
    pcomp <- pdat$pcomp
    up2 <- pdat$up2
    description <- pdat$description
    dat <- NULL
  } else stop(paste("secreted dataset", dataset, "not available"))
  print(paste0("pdat_secreted: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20191120
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, basis=basis, pcomp=pcomp, up2=up2, description=description))
}
