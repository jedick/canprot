# canprot/R/pdat_3D.R
# retrieve protein IDs for 3D cell culture (including tumor spheroid) datasets
# 20191125 initial version: some datasets moved from pdat_hypoxia.R
# 20191125-20191226 add data for 2020 compilation

pdat_3D <- function(dataset = 2020, basis = "rQEC") {
  if(identical(dataset, 2020)) {
    return(c("PLC+10",
             "MHG+12_P5", "MHG+12_P2",
             "MVC+12_perinecrotic", "MVC+12_necrotic",
             "YYW+13",
             "HKX+14",
             "RKP+14", "WRK+14",
             "MTK+15",
             "SSPR16=transcriptome", "YLW+16",
             "PPM+17=transcriptome",
             "KJK+18", "TGD18_NHF", "TGD18_CAF",
             "GSL+19", "HLC19", "LPK+19_preadipocytes", "LPK+19_adipocytes", "LPK+19_macrophages"
             ))
  }
  # remove tags
  dataset <- strsplit(dataset, "=")[[1]][1]
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse = "_")
  extdatadir <- system.file("extdata", package = "canprot")
  datadir <- paste0(extdatadir, "/expression/3D/")
  if(study=="MVC+12") {
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
    pcomp <- protcomp(dat$Entry, basis)
  } else if(study=="MHG+12") {
    # 20160415 MCF-7 tumourspheres, Morrison et al., 2012
    # MHG+12_P5, MHG+12_P2
    dat <- read.csv(paste0(datadir, "MHG+12.csv.xz"), as.is = TRUE)
    description <- paste("MCF-7", stage)
    # use data for specified experiment
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol[1]]), ]
    dat <- check_IDs(dat, "Protein.IDs")
    up2 <- dat[, icol[1]] < 0
    dat <- cleanup(dat, "Protein.IDs", up2)
    pcomp <- protcomp(dat$Protein.IDs, basis)
  } else if(study=="RKP+14") {
    # 20160718 organotypic spheroids, Rajcevic et al., 2014
    dat <- read.csv(paste0(datadir, "RKP+14.csv.xz"), as.is = TRUE)
    description <- "CRC-derived SPH"
    dat <- check_IDs(dat, "UniProt.Accession")
    up2 <- dat$Overall.Fold.Change > 0
    dat <- cleanup(dat, "UniProt.Accession", up2)
    pcomp <- protcomp(dat$UniProt.Accession, basis)
  } else if(study=="WRK+14") {
    # 20160721 3D spheroids / 2D culture, Wrzesinski et al., 2014
    dat <- read.csv(paste0(datadir, "WRK+14.csv.xz"), as.is = TRUE)
    description <- "HepG2/C3A SPH"
    # select highly changed proteins
    dat <- dat[!is.na(dat$log2.fold.change), ]
    dat <- dat[abs(dat$log2.fold.change) > 1, ]
    dat <- check_IDs(dat, "Entry")
    pcomp <- protcomp(dat$Entry, basis)
    up2 <- dat$log2.fold.change > 0
  } else if(study=="YLW+16") {
    # 20161109 HT29 colon cancer cell 3D/2D, Yue et al., 2011
    dat <- read.csv(paste0(datadir, "YLW+16.csv.xz"), as.is = TRUE)
    description <- "HT29 SPH"
    # find known UniProt IDs
    dat <- check_IDs(dat, "ProteinID")
    pcomp <- protcomp(dat$ProteinID, basis)
    up2 <- dat$Log2Rep1 > 0
  } else if(study=="TGD18") {
    # 20191125 cancer-associated fibroblasts, Tölle et al., 2018
    # TGD18_NHF, TGD18_CAF
    dat <- read.csv(paste0(datadir, "TGD18.csv.xz"), as.is = TRUE)
    if(stage == "NHF") description <- "normal human skin fibroblasts"
    if(stage == "CAF") description <- "cancer-associated fibroblasts"
    dat <- check_IDs(dat, "Protein.IDs")
    icol <- grep(stage, colnames(dat))
    up2 <- dat[, icol] > 0
    dat <- cleanup(dat, "Protein.IDs", up2)
    pcomp <- protcomp(dat$Protein.IDs, basis)
  } else if(study=="KJK+18") {
    # 20191126 SW480 cells, Kim et al., 2018
    dat <- read.csv(paste0(datadir, "KJK+18.csv.xz"), as.is = TRUE)
    description <- "SW480 cells"
    dat <- check_IDs(dat, "Majority.protein.IDs")
    pcomp <- protcomp(dat$Majority.protein.IDs, basis)
    up2 <- dat$Log2.3D.culture.2D.culture > 0
  } else if(study=="GSL+19") {
    # 20191127 4T1 cells (dense/sparse), Guo et al., 2019
    dat <- read.csv(paste0(datadir, "GSL+19.csv.xz"), as.is = TRUE)
    description <- "4T1 cells (dense/sparse)"
    dat <- check_IDs(dat, "Protein.ID", aa_file = paste0(extdatadir, "/aa/mouse/GSL+19_aa.csv"))
    pcomp <- protcomp(dat$Protein.ID, basis, aa_file = paste0(extdatadir, "/aa/mouse/GSL+19_aa.csv"))
    up2 <- dat$Ratio..Sparse.Dense. < 1
  } else if(study=="HKX+14") {
    # 20191206 U251 cells, He et al., 2014
    dat <- read.csv(paste0(datadir, "HKX+14.csv.xz"), as.is = TRUE)
    description <- "U251 cells"
    dat <- check_IDs(dat, "UniprotKB.AC")
    up2 <- dat$Ratio > 1
    pcomp <- protcomp(dat$UniprotKB.AC, basis)
  } else if(study=="HLC19") {
    # 20191206 HepG2 cells, Hurrell et al., 2019
    dat <- read.csv(paste0(datadir, "HLC19.csv.xz"), as.is = TRUE)
    description <- "HepG2 cells"
    up2 <- dat$Difference > 0
    pcomp <- protcomp(dat$Main.Accession, basis)
  } else if(study=="LPK+19") {
    # 20191206 3T3-L1 cells, Lee et al., 2019
    # LPK+19_preadipocytes, LPK+19_adipocytes, LPK+19_macrophages
    dat <- read.csv(paste0(datadir, "LPK+19.csv.xz"), as.is = TRUE)
    description <- paste("3T3-L1", stage)
    icol <- match(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    dat <- check_IDs(dat, "Accession", aa_file = paste0(extdatadir, "/aa/mouse/LPK+19_aa.csv"))
    up2 <- dat[, icol] > 0
    pcomp <- protcomp(dat$Accession, basis, aa_file = paste0(extdatadir, "/aa/mouse/LPK+19_aa.csv"))
  } else if(study=="MTK+15") {
    # 20191207 OV-90AD multicellular aggregates, Musrap et al., 2015
    dat <- read.csv(paste0(datadir, "MTK+15.csv.xz"), as.is = TRUE)
    description <- "OV-90AD multicellular aggregates"
    dat <- check_IDs(dat, "UniProt")
    up2 <- dat$Ratio.H.L.normalized > 1
    pcomp <- protcomp(dat$UniProt, basis)
  } else if(study=="PLC+10") {
    # 20191207 HepG2 cells, Pruksakorn et al., 2010
    dat <- read.csv(paste0(datadir, "PLC+10.csv.xz"), as.is = TRUE)
    description <- "HepG2 cells"
    dat <- check_IDs(dat, "Accession.no.")
    up2 <- dat$Difference == "Up"
    pcomp <- protcomp(dat$Accession.no., basis)
  } else if(study=="PPM+17") {
    # 20191207 HEY cell transcriptome, Paullin et al., 2017
    dat <- read.csv(paste0(datadir, "PPM+17.csv.xz"), as.is = TRUE)
    description <- "HEY cells transcriptome"
    dat <- check_IDs(dat, "Entry")
    up2 <- dat$Fold.change > 0
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis)
  } else if(study=="YYW+13") {
    # 20191207 HNSCC tumor spheres, Yan et al., 2013
    dat <- read.csv(paste0(datadir, "YYW+13.csv.xz"), as.is = TRUE)
    description <- "HNSCC tumor spheres"
    up2 <- dat$fold.change > 1
    pcomp <- protcomp(dat$Protein.ID, basis)
  } else if(study=="SSPR16") {
    # 20191226 HNSCC spheroids gene expression, Schmidt et al., 2016
    dat <- read.csv(paste0(datadir, "SSPR16.csv.xz"), as.is = TRUE)
    description <- "HNSCC spheroids transcriptome"
    up2 <- dat$logFC > 0
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis)
  } else stop(paste("3D dataset", dataset, "not available"))
  print(paste0("pdat_3D: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20190407
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset = dataset, basis = basis, pcomp = pcomp, up2 = up2, description = description))
}
