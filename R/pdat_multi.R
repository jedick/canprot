# canprot/R/pdat_multi.R
# retrieve protein IDs from studies that have multiple conditions
# 20191204 extracted from pdat_secreted.R

pdat_multi <- function(dataset = 2020, basis = "rQEC") {
  if(identical(dataset, 2020)) {
    return(c("CGH+17_exosomes=SEC", "CGH+17_secretome=SEC", "CGH+17_whole",
             "CLY+18_proteome", "CLY+18_secretome=SEC"
             ))
  }
  # remove tags
  dataset <- strsplit(dataset, "=")[[1]][1]
  # get study and stage/condition
  study <- strsplit(dataset, "_")[[1]][1]
  stage <- paste(strsplit(dataset, "_")[[1]][-1], collapse="_")
  extdatadir <- system.file("extdata", package="canprot")
  datadir <- paste0(extdatadir, "/expression/multi/")
  if(study=="CGH+17") {
    # 20190324 mouse cardiac fibroblast exosomes, secretome, whole-cell lysate, Cosme et al., 2017
    # CGH+17_exosomes, CGH+17_secretome, CGH+17_whole
    dat <- read.csv(paste0(datadir, "CGH+17.csv.xz"), as.is=TRUE)
    description <- paste("mouse CF", stage)
    # use selected dataset
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol[1]]), ]
    pcomp <- protcomp(dat$Entry, basis=basis, aa_file=paste0(extdatadir, "/aa/mouse/CGH+17_aa.csv"))
    up2 <- dat[, icol[1]] > 0
  } else if(study=="CLY+18") {
    # 20190324 HCT116 cells, Chen et al., 2018
    # CLY+18_proteome, CLY+18_secretome
    dat <- read.csv(paste0(datadir, "CLY+18.csv.xz"), as.is=TRUE)
    description <- paste("HCT116", stage)
    # use selected dataset
    icol <- grep(stage, colnames(dat))
    dat <- dat[!is.na(dat[, icol]), ]
    pcomp <- protcomp(dat$Accession, basis=basis)
    up2 <- dat[, icol] > 0
  } else stop(paste("multi dataset", dataset, "not available"))
  print(paste0("pdat_multi: ", description, " [", dataset, "]"))
  # use the up2 from the cleaned-up data, if it exists 20191120
  if("up2" %in% colnames(dat)) up2 <- dat$up2
  return(list(dataset=dataset, basis=basis, pcomp=pcomp, up2=up2, description=description))
}
