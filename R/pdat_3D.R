# canprot/R/pdat_3D.R
# retrieve protein IDs for 3D cell culture (including tumor spheroid) datasets
# 20191125 initial version: some datasets moved from pdat_hypoxia.R
# 20191125-20191230 add data for 2020 compilation

pdat_3D <- function(dataset = 2020, basis = "rQEC") {
  if(identical(dataset, 2020)) {
    return(c("PLC+10=cancer",
             "MHG+12_P5=cancer", "MHG+12_P2=cancer",
             "MVC+12_perinecrotic=cancer", "MVC+12_necrotic=cancer",
             "YYW+13=cancer", "ZMH+13_Matr.12h", "ZMH+13_Matr.24h",
             "HKX+14=cancer",
             "RKP+14=cancer", "WRK+14=cancer",
             "MTK+15=cancer",
             "MDN+16=transcriptome", "SSPR16=transcriptome=cancer", "YLW+16=cancer",
             "PPM+17=transcriptome=cancer",
             "KJK+18=cancer", "TGD18_NHF", "TGD18_CAF=cancer",
             "GADS19",
             "HLC19=cancer", "LPK+19_preadipocytes", "LPK+19_adipocytes", "LPK+19_macrophages"
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
    description <- paste("HT29 colon cancer cells", stage)
    if(stage=="perinecrotic") {
      # select proteins significantly changed in the perinecrotic region
      iPN <- dat$median.116.114 < 0.77 | dat$median.116.114 > 1.3
      dat <- dat[iPN, ]
      up2 <- dat$median.116.114 > 1.3
    }
    if(stage=="necrotic") {
      # select proteins significantly changed in the necrotic core
      iPN <- dat$median.117.114 < 0.77 | dat$median.117.114 > 1.3
      dat <- dat[iPN, ]
      up2 <- dat$median.117.114 > 1.3
    }
    pcomp <- protcomp(dat$Entry, basis)
  } else if(study=="MHG+12") {
    # 20160415 MCF-7 tumourspheres, Morrison et al., 2012
    # MHG+12_P5, MHG+12_P2
    dat <- read.csv(paste0(datadir, "MHG+12.csv.xz"), as.is = TRUE)
    description <- paste("MCF-7 breast cancer cells", stage)
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
    description <- "colorectal cancer-derived cells"
    dat <- check_IDs(dat, "UniProt.Accession")
    up2 <- dat$Overall.Fold.Change > 0
    dat <- cleanup(dat, "UniProt.Accession", up2)
    pcomp <- protcomp(dat$UniProt.Accession, basis)
  } else if(study=="WRK+14") {
    # 20160721 3D spheroids / 2D culture, Wrzesinski et al., 2014
    dat <- read.csv(paste0(datadir, "WRK+14.csv.xz"), as.is = TRUE)
    description <- "HepG2/C3A hepatocellular carcinoma"
    # select highly changed proteins
    dat <- dat[!is.na(dat$log2.fold.change), ]
    dat <- dat[abs(dat$log2.fold.change) > 1, ]
    dat <- check_IDs(dat, "Entry")
    pcomp <- protcomp(dat$Entry, basis)
    up2 <- dat$log2.fold.change > 0
  } else if(study=="YLW+16") {
    # 20161109 HT29 colon cancer cell 3D/2D, Yue et al., 2011
    dat <- read.csv(paste0(datadir, "YLW+16.csv.xz"), as.is = TRUE)
    description <- "HT29 colon carcinoma"
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
    description <- "SW480 colorectal cancer"
    dat <- check_IDs(dat, "Majority.protein.IDs")
    pcomp <- protcomp(dat$Majority.protein.IDs, basis)
    up2 <- dat$Log2.3D.culture.2D.culture > 0
  } else if(study=="HKX+14") {
    # 20191206 U251 cells, He et al., 2014
    dat <- read.csv(paste0(datadir, "HKX+14.csv.xz"), as.is = TRUE)
    description <- "U251 glioma cells"
    dat <- check_IDs(dat, "UniprotKB.AC")
    up2 <- dat$Ratio > 1
    pcomp <- protcomp(dat$UniprotKB.AC, basis)
  } else if(study=="HLC19") {
    # 20191206 HepG2 cells, Hurrell et al., 2019
    dat <- read.csv(paste0(datadir, "HLC19.csv.xz"), as.is = TRUE)
    description <- "HepG2 hepatocellular carcinoma"
    # 20191230 this should be less than to make the differences consistent with text (590 up, 573 down)
    up2 <- dat$Difference < 0
    pcomp <- protcomp(dat$Main.Accession, basis)
  } else if(study=="LPK+19") {
    # 20191206 3T3-L1 cells, Lee et al., 2019
    # LPK+19_preadipocytes, LPK+19_adipocytes, LPK+19_macrophages
    dat <- read.csv(paste0(datadir, "LPK+19.csv.xz"), as.is = TRUE)
    description <- paste("3T3-L1", stage)
    icol <- match(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    dat <- check_IDs(dat, "Accession", aa_file = paste0(extdatadir, "/aa/mouse/LPK+19_aa.csv.xz"))
    up2 <- dat[, icol] > 0
    pcomp <- protcomp(dat$Accession, basis, aa_file = paste0(extdatadir, "/aa/mouse/LPK+19_aa.csv.xz"))
  } else if(study=="MTK+15") {
    # 20191207 OV-90AD multicellular aggregates, Musrap et al., 2015
    dat <- read.csv(paste0(datadir, "MTK+15.csv.xz"), as.is = TRUE)
    description <- "OV-90AD ovarian cancer multicellular aggregates"
    dat <- check_IDs(dat, "UniProt")
    up2 <- dat$Ratio.H.L.normalized > 1
    pcomp <- protcomp(dat$UniProt, basis)
  } else if(study=="PLC+10") {
    # 20191207 HepG2 cells, Pruksakorn et al., 2010
    dat <- read.csv(paste0(datadir, "PLC+10.csv.xz"), as.is = TRUE)
    description <- "HepG2 hepatocellular carcinoma cells"
    dat <- check_IDs(dat, "Accession.no.")
    up2 <- dat$Difference == "Up"
    pcomp <- protcomp(dat$Accession.no., basis)
  } else if(study=="PPM+17") {
    # 20191207 HEY cell transcriptome, Paullin et al., 2017
    dat <- read.csv(paste0(datadir, "PPM+17.csv.xz"), as.is = TRUE)
    description <- "HEY epithelial ovarian cancer transcriptome"
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
  } else if(study=="GADS19") {
    # 20191230 skin fibroblasts 3D/2D with various treatments, Gęgotek et al., 2019
    dat <- read.csv(paste0(datadir, "GADS19.csv.xz"), as.is = TRUE)
    description <- "skin fibroblasts 3D / 2D"
    dat <- check_IDs(dat, "ID")
    up2 <- rep(NA, nrow(dat))
    # proteins with differential expression in at least 4 of 6 treatments
    up2[rowSums(dat[, -1] > 1.2) > 3] <- TRUE
    up2[rowSums(dat[, -1] < 1/1.2) > 3] <- FALSE
    dat <- cleanup(dat, "ID", up2)
    pcomp <- protcomp(dat$ID, basis)
  } else if(study=="ZMH+13") {
    # 20200114 human umbilical vein endothelial cells, Zanivan et al., 2013
    # grown in Matrigel (3D): ZMH+13_Matr.12h, ZMH+13_Matr.24h, ZMH+13_Matr.GFR, ZMH+13_Matr.30h
    # spreading on culture dish: ZMH+13_Matr.dil, ZMH+13_LAM, ZMH+13_FN, ZMH+13_BSA
    dat <- read.csv(paste0(datadir, "ZMH+13.csv.xz"), as.is = TRUE)
    description <- paste("HUVEC", gsub("Matr.", "Matrigel ", stage, fixed = TRUE))
    # get data for selected experiment
    idn <- match(paste0(stage, ".down"), colnames(dat))
    iup <- match(paste0(stage, ".up"), colnames(dat))
    dat <- dat[dat[, idn] == "yes" | dat[, iup] == "yes", ]
    dat <- check_IDs(dat, "Protein.IDs")
    pcomp <- protcomp(dat$Protein.IDs, basis)
    up2 <- dat[, iup] == "yes"
  } else if(study=="MDN+16") {
    # 20200115 placental syncytiotrophoblast JEG-3 cells transcriptome, McConkey et al., 2016
    dat <- read.csv(paste0(datadir, "MDN+16.csv.xz"), as.is = TRUE)
    description <- "placental syncytiotrophoblast JEG-3 cells transcriptome"
    dat <- check_IDs(dat, "Entry")
    up2 <- dat$log2FoldChange > 0
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis)
  } else stop(paste("3D dataset", dataset, "not available"))
  print(paste0("pdat_3D: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20190407
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset = dataset, basis = basis, pcomp = pcomp, up2 = up2, description = description))
}
