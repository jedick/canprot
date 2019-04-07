# canprot/R/zzz.R
# create the data objects for amino acid compositions of human proteins
# 20160705 jmd (data/canprot.R)
# 20190225 moved to R/zzz.R

# the canprot environment is made here in open code 20190214
# https://stackoverflow.com/questions/41954302/where-to-create-package-environment-variables
canprot <- new.env()

# initialize the canprot environment
.onAttach <- function(libname,pkgname) {
  with(canprot, {
    # read amino acid compositions of human proteins and show some information
    #human_base <- read.csv("human_base.csv.xz", as.is=TRUE)
    load(system.file("/extdata/protein/human_base.Rdata", package="canprot"))
    packageStartupMessage(paste("human_base:", nrow(human_base), "proteins"))
    #human_additional <- read.csv("human_additional.csv.xz", as.is=TRUE)
    load(system.file("/extdata/protein/human_additional.Rdata", package="canprot"))
    packageStartupMessage(paste("human_additional:", nrow(human_additional), "proteins"))
    human_extra <- read.csv(system.file("/extdata/protein/human_extra.csv", package="canprot"), as.is=TRUE)
    packageStartupMessage(paste("human_extra:", nrow(human_extra), "proteins"))
    # create the object
    human_aa <- rbind(human_base, human_additional, human_extra)
    # warn if there are duplicated proteins
    local({
      idup <- duplicated(human_aa$protein)
      if(any(idup)) warning("data(human): duplicated proteins: ", paste(human_aa$protein[idup], collapse=" "))
    })
    # updates for UniProt IDs
    uniprot_updates <- read.csv(system.file("/extdata/protein/uniprot_updates.csv", package="canprot"), as.is=TRUE)
  })
}
