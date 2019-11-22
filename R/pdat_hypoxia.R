# canprot/R/pdat_hypoxia.R
# retrieve protein IDs for hypoxia data sets
# 20160414 jmd

pdat_hypoxia <- function(dataset=NULL, basis="rQEC") {
  if(is.null(dataset)) {
    return(c("HXS+06",
             "BRA+10", "DPL+10",
             "BMJ+11", "CBW+11",
             "LAR+12", "MHG+12_P5=SPH", "MHG+12_P2=SPH",
             "MVC+12_perinecrotic=SPH", "MVC+12_necrotic=SPH",
             "FWH+13", "RHD+13_Hx48", "RHD+13_Hx72", "RHD+13_ReOx", "VTMF13",
             "DYL+14_Hx48-S", "DYL+14_Hx72-S", "DYL+14_ReOx-S",
             "DYL+14_Hx48-P", "DYL+14_Hx72-P", "DYL+14_ReOx-P",
             "RKP+14=SPH", "WRK+14=SPH",
             "BSA+15",
             "HWA+16",
             "LCS16_transcription", "LCS16_translation",
             "RSE+16=ASC", "XCJ+16_CoCl2", "XCJ+16_SAL=ReOx", "YLW+16=SPH"))
  }
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse="_")
  extdatadir <- system.file("extdata", package="canprot")
  datadir <- paste0(extdatadir, "/expression/hypoxia/")
  if(study=="BSA+15") {
    # 20160412 HeLa hypoxia, Bousquet et al., 2005
    dat <- read.csv(paste0(datadir, "BSA+15.csv.xz"), as.is=TRUE)
    description <- "HeLa"
    # use updated UniProt IDs
    inew <- which(dat$Uniprot.new != "")
    dat$Uniprot.accession[inew] <- dat$Uniprot.new[inew]
    # remove Q8IWE2, which is duplicated with opposite expression ratio
    dat <- dat[dat$Uniprot.accession!="Q8IWE2", ]
    pcomp <- protcomp(dat$Uniprot.accession, basis=basis)
    up2 <- dat$Ratio..H.L. > 1
  } else if(study=="MVC+12") {
    # 20160413 spheriod hypoxia, McMahon et al., 2012
    # MVC+12_perinecrotic, MVC+12_necrotic
    dat <- read.csv(paste0(datadir, "MVC+12.csv.xz"), as.is=TRUE)
    if(stage=="perinecrotic") {
      # select proteins significantly changed in the perinecrotic region
      iPN <- dat$median.116.114 < 0.77 | dat$median.116.114 > 1.3
      dat <- dat[iPN, ]
      up2 <- dat$median.116.114 > 1.3
      description <- "SPH perinecrotic"
    }
    if(stage=="necrotic") {
      # select proteins significantly changed in the necrotic core
      iPN <- dat$median.117.114 < 0.77 | dat$median.117.114 > 1.3
      dat <- dat[iPN, ]
      up2 <- dat$median.117.114 > 1.3
      description <- "SPH necrotic"
    }
    pcomp <- protcomp(dat$Entry, basis=basis)
  } else if(study=="MHG+12") {
    # 20160415 MCF-7 tumourspheres, Morrison et al., 2012
    # MHG+12_P5, MHG+12_P2
    dat <- read.csv(paste0(datadir, "MHG+12.csv.xz"), as.is=TRUE)
    description <- paste("MCF-7 SPH", stage)
    # use dat for specified experiment
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol[1]]), ]
    dat <- check_IDs(dat, "Protein.IDs")
    up2 <- dat[, icol[1]] < 0
    dat <- cleanup(dat, "Protein.IDs", up2)
    pcomp <- protcomp(dat$Protein.IDs, basis=basis)
  } else if(study=="HXS+06") {
    # 20160415 leukemic U937 cells, Han et al., 2006
    dat <- read.csv(paste0(datadir, "HXS+06.csv.xz"), as.is=TRUE)
    description <- "U937"
    # update IDs with new ones
    inew <- dat$UniProt.new != ""
    dat$Swiss.Prot.accession.no.[inew] <- dat$UniProt.new[inew]
    up2 <- dat$Mean.fold.H.N > 1
    dat <- cleanup(dat, "Swiss.Prot.accession.no.", up2)
    pcomp <- protcomp(dat$Swiss.Prot.accession.no., basis=basis)
  } else if(study=="RHD+13") {
    # 20160419 A431 cells, Ren et al., 2013
    # RHD+13_Hx48, RHD+13_Hx72, RHD+13_ReOx
    dat <- read.csv(paste0(datadir, "RHD+13.csv.xz"), as.is=TRUE)
    description <- paste("A431", stage)
    # columns with the ratios and p-values
    if(stage=="Hx48") icol <- grep("115", colnames(dat))
    if(stage=="Hx72") icol <- grep("116", colnames(dat))
    if(stage=="ReOx") icol <- grep("117", colnames(dat))
    # use significantly highly changed proteins
    dat <- dat[dat[, icol[2]] < 0.05, ]
    dat <- dat[dat[, icol[1]] > sqrt(2) | dat[, icol[1]] < 1/sqrt(2), ]
    # get known UniProt IDs
    dat <- check_IDs(dat, "Accession")
    up2 <- dat[, icol[1]] > 1
    dat <- cleanup(dat, "Accession", up2)
    pcomp <- protcomp(dat$Accession, basis=basis)
  } else if(study=="BMJ+11") {
    # 20160713 DU145 cells prolonged hypoxia, van den Beucken et al., 2011
    dat <- read.csv(paste0(datadir, "BMJ+11.csv.xz"), as.is=TRUE)
    description <- "DU145"
    # keep proteins detected in prolonged hypoxia
    dat <- dat[dat$induced_prolonged | dat$repressed_prolonged, ]
    up2 <- dat$induced_prolonged
    dat <- cleanup(dat, "uniprot", up2)
    pcomp <- protcomp(dat$uniprot, basis=basis)
  } else if(study=="FWH+13") {
    # 20160716 THP-1 macrophages CV (control virus) hypoxia, Fuhrmann et al., 2013
    dat <- read.csv(paste0(datadir, "FWH+13.csv.xz"), as.is=TRUE)
    description <- "THP-1"
    up2 <- dat$Norm.CH > 0
    pcomp <- protcomp(dat$UniProt, basis=basis)
  } else if(study=="HWA+16") {
    # 20160716 U87MG and 786-O translatome, Ho et al., 2016
    dat <- read.csv(paste0(datadir, "HWA+16.csv.xz"), as.is=TRUE)
    # keep those with fold change < 0.5 or > 2 (log2 < -1 or > 1)
    dat <- dat[ dat$Hypoxia.Heavy - dat$Normoxia.Heavy > 1 |
                dat$Hypoxia.Heavy - dat$Normoxia.Heavy < -1, ]
    description <- "U87MG and 786-O"
    up2 <- dat$Hypoxia.Heavy - dat$Normoxia.Heavy > 0
    dat <- check_IDs(dat, "Uniprot.Accession")
    pcomp <- protcomp(dat$Uniprot.Accession, basis=basis)
  } else if(study=="RKP+14") {
    # 20160718 organotypic spheroids, Rajcevic et al., 2014
    dat <- read.csv(paste0(datadir, "RKP+14.csv.xz"), as.is=TRUE)
    description <- "CRC-derived SPH"
    dat <- check_IDs(dat, "UniProt.Accession")
    up2 <- dat$Overall.Fold.Change > 0
    dat <- cleanup(dat, "UniProt.Accession", up2)
    pcomp <- protcomp(dat$UniProt.Accession, basis=basis)
  } else if(study=="CBW+11") {
    # 20160720 neuroblastoma cells, Cifani et al., 2011
    dat <- read.csv(paste0(datadir, "CBW+11.csv.xz"), as.is=TRUE)
    description <- "SK-N-BE(2)c; IMR-32"
    # remove "-1" suffix for isoform 1
    dat$UniProt <- gsub("-1", "", dat$UniProt)
    dat <- check_IDs(dat, "UniProt")
    up2 <- dat$Be2c > 1
    dat <- cleanup(dat, "UniProt", up2)
    pcomp <- protcomp(dat$UniProt, basis=basis)
  } else if(study=="WRK+14") {
    # 20160721 3D spheroids / 2D culture, Wrzesinski et al., 2014
    dat <- read.csv(paste0(datadir, "WRK+14.csv.xz"), as.is=TRUE)
    description <- "HepG2/C3A SPH"
    # select highly changed proteins
    dat <- dat[!is.na(dat$log2.fold.change), ]
    dat <- dat[abs(dat$log2.fold.change) > 1, ]
    dat <- check_IDs(dat, "Entry")
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- dat$log2.fold.change > 0
  } else if(study=="DPL+10") {
    # 20160722 B104 rat neuroblastoma cells, Datta et al., 2010
    dat <- read.csv(paste0(datadir, "DPL+10.csv.xz"), as.is=TRUE)
    description <- "B104"
    # select highly changed proteins
    dat <- dat[dat$HYP.LSC > 1.2 | dat$HYP.LSC < 0.83, ]
    dat <- check_IDs(dat, "UniProt", aa_file=paste0(extdatadir, "/aa/rat/DPL+10_aa.csv"))
    up2 <- dat$HYP.LSC > 1
    dat <- cleanup(dat, "UniProt", up2)
    pcomp <- protcomp(dat$UniProt, basis=basis, aa_file=paste0(extdatadir, "/aa/rat/DPL+10_aa.csv"))
  } else if(study=="LCS16") {
    # 20160728 HCT116 transcription and translation, Lai et al., 2016
    # LCS16_transcription, LCS16_translation
    dat <- read.csv(paste0(datadir, "LCS16.csv.xz"), as.is=TRUE)
    description <- paste("HCT116", stage)
    # select the experiment
    icol <- grep(stage, tolower(colnames(dat)))
    idiff <- sapply(dat[, icol[1]] | dat[, icol[2]], isTRUE)
    dat <- dat[idiff, ]
    # which proteins (genes) are up-regulated
    iicol <- icol[grep("up", colnames(dat)[icol])]
    up2 <- sapply(dat[, iicol], isTRUE)
    dat <- cleanup(dat, "UniProt", up2)
    pcomp <- protcomp(dat$UniProt, basis=basis)
  } else if(study=="DYL+14") {
    # 20160729 A431 cells, Dutta et al., 2014
    # DYL+14_Hx48-S, DYL+14_Hx72-S, DYL+14_ReOx-S,
    # DYL+14_Hx48-P, DYL+14_Hx72-P, DYL+14_ReOx-P
    dat <- read.csv(paste0(datadir, "DYL+14.csv.xz"), as.is=TRUE)
    description <- paste("A431", stage)
    # -S (supernatant) and -P (pellet) datasets
    if(stage=="Hx48-S") icol <- grep("114", colnames(dat))
    if(stage=="Hx72-S") icol <- grep("115", colnames(dat))
    if(stage=="ReOx-S") icol <- grep("116", colnames(dat))
    if(stage=="Hx48-P") icol <- grep("118", colnames(dat))
    if(stage=="Hx72-P") icol <- grep("119", colnames(dat))
    if(stage=="ReOx-P") icol <- grep("121", colnames(dat))
    # now deal with -P (pellet) datasets
    # divide (expt / Nx-S) by (Nx-P / Nx-S)
    if(grepl("-P", stage)) dat[, icol[1]] <- dat[, icol[1]] / dat$X117.113
    # keep significantly changed proteins based on p-value and ratio
    if(!grepl("-P", stage)) dat <- dat[dat[, icol[2]] < 0.05, ]
    dat <- dat[dat[, icol[1]] > sqrt(2) | dat[, icol[1]] < 1/sqrt(2), ]
    dat <- check_IDs(dat, "Accession")
    pcomp <- protcomp(dat$Accession, basis=basis)
    up2 <- dat[, icol[1]] > 1
  } else if(study=="RSE+16") {
    # 20160729 adipose-derived stem cells, Riis et al., 2016
    dat <- read.csv(paste0(datadir, "RSE+16.csv.xz"), as.is=TRUE)
    description <- "adipose-derived SC"
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- dat$Regulated == "up"
  } else if(study=="VTMF13") {
    # 20160804 neuroblastoma cell line, Villeneuve et al., 2013
    dat <- read.csv(paste0(datadir, "VTMF13.csv.xz"), as.is=TRUE)
    description <- "SH-SY5Y"
    # keep proteins with large expression ratio
    dat <- dat[dat$Ratio.H.L.Normalized > 1.2 | dat$Ratio.H.L.Normalized < 0.83, ]
    # find known UniProt IDs
    dat <- check_IDs(dat, "Uniprot")
    up2 <- dat$Ratio.H.L.Normalized > 1.2
    dat <- cleanup(dat, "Uniprot", up2)
    pcomp <- protcomp(dat$Uniprot, basis=basis)
  } else if(study=="BRA+10") {
    # 20160805 placental tissue secretome, Blankley et al., 2010
    dat <- read.csv(paste0(datadir, "BRA+10.csv.xz"), as.is=TRUE)
    description <- "placental secretome"
    dat$UniProt.accession <- sapply(strsplit(dat$UniProt.accession, "|", fixed=TRUE), "[", 2)
    dat <- check_IDs(dat, "UniProt.accession")
    pcomp <- protcomp(dat$UniProt.accession, basis=basis)
    up2 <- dat$Fold.change > 0
  } else if(study=="LAR+12") {
    # 20160826 rat heart ischemia, Li et al., 2012
    dat <- read.csv(paste0(datadir, "LAR+12.csv.xz"), as.is=TRUE)
    description <- "H9C2"
    up2 <- dat$Isch.Ctrl > 1
    dat <- cleanup(dat, "UniProt", up2)
    pcomp <- protcomp(dat$UniProt, basis=basis, aa_file=paste0(extdatadir, "/aa/mouse/LAR+12_aa.csv"))
  } else if(study=="YLW+16") {
    # 20161109 HT29 colon cancer cell 3D/2D, Yue et al., 2011
    dat <- read.csv(paste0(datadir, "YLW+16.csv.xz"), as.is=TRUE)
    description <- "HT29 SPH"
    # find known UniProt IDs
    dat <- check_IDs(dat, "ProteinID")
    pcomp <- protcomp(dat$ProteinID, basis=basis)
    up2 <- dat$Log2Rep1 > 0
  } else if(study=="XCJ+16") {
    # 20161119 cardiomyocytes CoCl2 (hypoxia mimetic) or SAL (anti-hypoxic), Xu et al., 2016
    # XCJ+16_CoCl2, XCJ+16_SAL
    dat <- read.csv(paste0(datadir, "XCJ+16.csv.xz"), as.is=TRUE)
    description <- paste("cardiomyocytes", stage)
    # use selected dataset
    icol <- grep(paste0("Log2.", stage), colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    pcomp <- protcomp(dat$Entry, basis=basis, aa_file=paste0(extdatadir, "/aa/rat/XCJ+16_aa.csv"))
    up2 <- dat[, icol] > 0
  } else stop(paste("hypoxia dataset", dataset, "not available"))
  print(paste0("pdat_hypoxia: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20191120
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, basis=basis, pcomp=pcomp, up2=up2, description=description))
}
