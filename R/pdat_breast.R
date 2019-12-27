# canprot/R/pdat_breast.R
# retrieve protein IDs for breast cancer studies
# 20160411-20191120 assemble data for 2020 compilation

pdat_breast <- function(dataset = 2020, basis = "rQEC") {
  if(identical(dataset, 2020)) {
    return(c(
             "AMG+08",
             "CIR+10=luminalA",
             "SRG+10",
             "HTP+11",
             "GTM+12_IDC.benign", "GTM+12_IDC",
             "LLL+13=basal", "SRS+13_DCIS", "SRS+13_IC",
             "GSB+14_tumor",
             "PPH+14",
             "CVJ+15=basal",
             "FKZ+15_A=luminalA=transcriptome", "FKZ+15_B=luminalB=transcriptome", "FKZ+15_TN=basal=transcriptome",
             "PGT+16_ES", "PGT+16_LS",
             "PBR+16_tumor",
             "BST+17_epithelium",
             "TZD+18_all", "TZD+18_basal",
             "GCS+19_PTxNCT", "GCS+19_PTxANT", "LLC+19"
             ))
  }
  # remove tags
  dataset <- strsplit(dataset, "=")[[1]][1]
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse="_")
  extdatadir <- system.file("extdata", package="canprot")
  datadir <- paste0(extdatadir, "/expression/breast/")
  if(study=="GTM+12") {
    # 20160411 higher in benign breast growth or IDC than normal adjacent tissue, Gormley et al., 2012
    # GTM+12_benign, GTM+12_IDC, GTM+12_IDC.benign
    dat <- read.csv(paste0(datadir, "GTM+12.csv.xz"), as.is=TRUE)
    if(stage=="benign") {
      description <- "benign / normal adj."
      dat <- dat[which(dat$Benign.vs..NormalAdjacent), ]
      up2 <- dat$Benign.Mean > dat$NormalAdjacent.Mean
    }
    if(stage=="IDC") {
      description <- "IDC / normal adj."
      dat <- dat[which(dat$NormalAdjacent.vs..Infiltrating.Ductal.Carcinoma), ]
      up2 <- dat$IDC.Mean > dat$NormalAdjacent.Mean
    }
    if(stage=="IDC.benign") {
      description <- "IDC / benign"
      dat <- dat[which(dat$Benign.vs..Infiltrating.Ductal.Carcinoma), ]
      up2 <- dat$IDC.Mean > dat$Benign.Mean
    }
    pcomp <- protcomp(dat$UniprotID, basis=basis)
  } else if(study=="CIR+10") {
    # 20160414 breast cancer Cha et al., 2010
    dat <- read.csv(paste0(datadir, "CIR+10.csv.xz"), as.is=TRUE)
    description <- "ER+ T/N"
    pcomp <- protcomp(dat$Uniprot, basis=basis)
    up2 <- dat[, "SpI"] > 0
  } else if(study=="PBR+16") {
    # 20160415 breast cancer Pozniak et al., 2016
    # PBR+16_tumor, PBR+16_LNP
    dat <- read.csv(paste0(datadir, "PBR+16.csv.xz"), as.is=TRUE)
    description <- "T/N"
    dat <- check_IDs(dat, "Protein.IDs")
    pcomp <- protcomp(dat$Protein.IDs, basis=basis)
    up2 <- dat$Cluster == "Upregulated in Tumor"
  } else if(study=="SRG+10") {
    # 20160415 breast cancer Sutton et al., 2010
    dat <- read.csv(paste0(datadir, "SRG+10.csv.xz"), as.is=TRUE)
    description <- "invasive carcinoma / normal"
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- dat$change == "increased"
  } else if(study=="HTP+11") {
    # 20160420 microvesicles Hill et al., 2011
    dat <- read.csv(paste0(datadir, "HTP+11.csv.xz"), as.is=TRUE)
    description <- "microvessels IDC / nonmalignant"
    # update IDs with new ones
    inew <- dat$UniProt.new != ""
    dat$Swiss.Prot[inew] <- dat$UniProt.new[inew]
    pcomp <- protcomp(dat$Swiss.Prot, basis=basis)
    up2 <- dat$Expressed == "over"
  } else if(study=="LLL+13") {
    # 20160420 triple negative breast cancer Liang et al., 2013
    dat <- read.csv(paste0(datadir, "LLL+13.csv.xz"), as.is=TRUE)
    description <- "TNBC tumor / paraneoplastic"
    up2 <- dat$Style == "up"
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis=basis)
  } else if(study=="PPH+14") {
    # 20160421 tumor subtypes vs normal, Panis et al., 2014
    dat <- read.csv(paste0(datadir, "PPH+14.csv.xz"), as.is=TRUE)
    description <- "multiple subtypes T/N"
    # keep proteins that are up or down in all subtypes
    dat <- dat[sapply(apply(dat[, 4:7], 1, unique), length)==1, ]
    up2 <- dat$TN == "up"
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis=basis)
  } else if(study=="SRS+13") {
    # 20160421 tumor vs normal, Shaheed et al., 2013
    # SRS+13_DCIS, SRS+13_IC
    # extract differential proteins from datafile 20170831
    dat <- read.csv(paste0(datadir, "SRS+13.csv.xz"), as.is=TRUE)
    description <- paste(stage, "/ normal")
    # remove proteins located to blood
    dat <- dat[!dat$Location=="blood", ]
    # remove ratio data with less than 2 peptides
    iratio <- grep("_ratio", colnames(dat))
    ipeptides <- grep("peptides", colnames(dat))
    dat[, ipeptides][is.na(dat[, ipeptides])] <- 0
    dat[, iratio][(dat[, ipeptides]) < 2] <- NA
    # normalize tumor / normal means for each patient
    # get the log ratios and their means
    logratio <- log(dat[,iratio])
    meanlog <- apply(logratio, 2, mean, na.rm=TRUE)
    # normalize the log values (subtract the mean difference)
    normlog <- t(t(logratio) - meanlog)
    # use the normalized log ratios
    dat[, iratio] <- round(10^normlog, 2)
    # where to save the the up/down classification
    dat <- cbind(dat, expression=NA)
    if(stage=="DCIS") {
      # get DCIS differences
      iDCIS <- match(c("pt3A_ratio", "pt7B_ratio", "pt9C_ratio"), colnames(dat))
      nup <- apply(dat[, iDCIS] > 1, 1, sum, na.rm=TRUE)
      ndown <- apply(dat[, iDCIS] < 1, 1, sum, na.rm=TRUE)
      iup <- nup > 1 & ndown == 0
      idown <- ndown > 1 & nup == 0
      dat$expression[iup] <- "up"
      dat$expression[idown] <- "down"
    } else if(stage=="IC") {
      # get IC differences
      iIC <- match(c("pt5A_ratio", "pt4A_ratio", "pt6B_ratio", "pt2B_ratio"), colnames(dat))
      nup <- apply(dat[, iIC] > 1, 1, sum, na.rm=TRUE)
      ndown <- apply(dat[, iIC] < 1, 1, sum, na.rm=TRUE)
      iup <- nup > 2 & ndown == 0
      idown <- ndown > 2 & nup == 0
      dat$expression[iup] <- "up"
      dat$expression[idown] <- "down"
    }
    # use proteins classified as up or down
    dat <- dat[!is.na(dat$expression), ]
    pcomp <- protcomp(dat$Entry, basis=basis)
    up2 <- dat$expression == "up"
  } else if(study=="PGT+16") {
    # 20160718 luminal B or HER enriched; early or late vs control, Pendharkar et al., 2016
    # PGT+16_LB, PGT+16_ER, PGT+16_LB, PGT+16_HE
    dat <- read.csv(paste0(datadir, "PGT+16.csv.xz"), as.is=TRUE)
    if(stage=="LB") description <- "Luminal B HER2 positive T / N"
    if(stage=="HE") description <- "HER2 enriched T / N"
    if(stage=="ES") description <- "early stage T / N"
    if(stage=="LS") description <- "late stage T / N"
    # just use proteins quantified at this stage
    icol <- grep(paste0("Fold.change.", stage), colnames(dat))
    dat <- check_IDs(dat, "Accession.number")
    up2 <- dat[, icol] > 0
    dat <- cleanup(dat, "Accession.number", up2)
    pcomp <- protcomp(dat$Accession.number, basis=basis)
  } else if(study=="GSB+14") {
    # 20170116 breast adenocarcinoma and treatment of ZR-75-1 with TGFβ or IL-1β, Groessl et al., 2014
    # GSB+14_tumor, GSB+14_near, GSB+14_IL-1β, GSB+14_TGFβ
    dat <- read.csv(paste0(datadir, "GSB+14.csv.xz"), as.is=TRUE, check.names = FALSE)
    if(stage %in% c("tumor", "near")) description <- paste0("adenocarcinoma ", stage, "/distant")
    #if(stage %in% c("IL-1β", "TGFβ")) description <- paste("breast ZR-75-1;", stage)
    # use specified experiment
    icol <- grep(stage, colnames(dat))
    if(stage %in% c("tumor", "near")) icon <- grep("LFQ distant", colnames(dat))
    #if(stage %in% c("IL-1β", "TGFβ")) icon <- grep("LFQ con", colnames(dat))
    # fold change > 1.05 or < 0.95 
    dat <- dat[!is.na(dat[, icol]), ]
    rat <- dat[, icol] / dat[, icon]
    dat <- dat[rat > 1.05 | rat < 0.95, ]
    dat <- check_IDs(dat, "Accession")
    up2 <- dat[, icol] / dat[, icon] > 1
    dat <- cleanup(dat, "Accession", up2)
    pcomp <- protcomp(dat$Accession, basis=basis)
  } else if(study=="BST+17") {
    # 20170811 breast epithelium and stroma, Braakman et al., 2017
    # BST+17_epithelium, BST+17_stroma
    dat <- read.csv(paste0(datadir, "BST+17.csv.xz"), as.is=TRUE)
    description <- paste("tumor", stage, "/ normal")
    # use selected tissue
    icol <- grep(stage, colnames(dat))
    # use significantly differentially abundant proteins
    dat <- dat[dat[, icol[2]] < 0.05, ]
    dat <- check_IDs(dat, "Protein.IDs")
    pcomp <- protcomp(dat$Protein.IDs, basis=basis)
    up2 <- dat[, icol[1]] > 0
  } else if(study=="CVJ+15") {
    # 20170814 TNBC tumor / normal, Campone et al., 2015
    dat <- read.csv(paste0(datadir, "CVJ+15.csv.xz"), as.is=TRUE)
    description <- "TNBC tumor / normal"
    dat <- check_IDs(dat, "Accession.Number")
    up2 <- dat$Mean > 1
    dat <- cleanup(dat, "Accession.Number", up2)
    pcomp <- protcomp(dat$Accession.Number, basis=basis)
  } else if(study=="AMG+08") {
    # 20170828 cancer / periphery, Alldridge et al., 2008
    dat <- read.csv(paste0(datadir, "AMG+08.csv.xz"), as.is=TRUE)
    description <- "cancer / periphery"
    # remove "-1" isoform suffixes
    dat$ID <- gsub("-1", "", dat$ID)
    dat <- check_IDs(dat, "ID")
    up2 <- dat$Cancer_Peptides > 0
    dat <- cleanup(dat, "ID", up2)
    pcomp <- protcomp(dat$ID, basis=basis)
  } else if(study=="FKZ+15") {
    # 20170829 gene expression luminal A, luminal B, triple negative / normal, Fu et al., 2015
    # FKZ+15_A, FKZ+15_B, FKZ+15_TN
    dat <- read.csv(paste0(datadir, "FKZ+15.csv.xz"), as.is=TRUE)
    if(stage=="A") { description <- "luminal A / normal transcriptome"; icol <- grep("luminalA", colnames(dat))[2] }
    if(stage=="B") { description <- "luminal B / normal transcriptome"; icol <- grep("luminalB", colnames(dat))[2] }
    if(stage=="TN") { description <- "triple negative / normal transcriptome"; icol <- grep("TN", colnames(dat))[2] }
    # use indicated subtype
    dat <- dat[!is.na(dat[, icol]), ]
    # remove "-1" isoform suffixes
    dat$Entry <- gsub("-1$", "", dat$Entry)
    dat <- check_IDs(dat, "Entry")
    up2 <- dat[, icol] > 0
    dat <- cleanup(dat, "Entry", up2)
    # isoform suffixes (other than -1) listed here in uniprot_updates correspond to canonical sequence (so use plain IDs)
    pcomp <- protcomp(dat$Entry, basis=basis)
  } else if(study=="LLC+19") {
    # 20190318 tumor / adjacent normal, Liu et al., 2019
    dat <- read.csv(paste0(datadir, "LLC+19.csv.xz"), as.is=TRUE)
    description <- "T / N"
    pcomp <- protcomp(dat$Protein.accession, basis=basis)
    up2 <- dat$Regulated.Type == "Up"
  } else if(study=="TZD+18") {
    # 20190321 tumor / adjacent normal (all subtypes or basal), Tang et al., 2018
    # TZD+18_all, TZD+18_basal
    dat <- read.csv(paste0(datadir, "TZD+18.csv.xz"), as.is=TRUE)
    description <- paste("T / N", stage)
    # keep proteins differentially expressed in all or basal tumors
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    dat <- check_IDs(dat, "UNIPROT")
    up2 <- dat[, icol] > 0
    dat <- cleanup(dat, "UNIPROT", up2)
    pcomp <- protcomp(dat$UNIPROT, basis=basis)
  } else if(study=="GCS+19") {
    # 20191120 tumor and lymph node vs contralateral and adjacent, Gomig et al., 2019
    # GCS+19_PTxNCT, GCS+19_LNxNCT, GCS+19_PTxLN, GCS+19_PTxANT, GCS+19_LNxANT
    dat <- read.csv(paste0(datadir, "GCS+19.csv.xz"), as.is=TRUE)
    description <- stage
    icol <- grep(stage, colnames(dat))
    up2 <- dat[, icol] > 0
    dat <- cleanup(dat, "Entry", up2)
    pcomp <- protcomp(dat$Entry, basis=basis)
  } else stop(paste("breast dataset", dataset, "not available"))
  print(paste0("pdat_breast: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20190429
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, pcomp=pcomp, up2=up2, description=description))
}
