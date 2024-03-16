# canprot/zzz.R
# Create the data objects for amino acid compositions of human proteins
# 20160705 first version (data/canprot.R)
# 20190225 moved to R/zzz.R
# 20200509 change environment name from canprot to human
# 20240301 revert environment name to canprot

# The 'canprot' environment is made here in open code 20190214
# https://stackoverflow.com/questions/41954302/where-to-create-package-environment-variables
canprot <- new.env()

# Initialize the 'canprot' environment
.onAttach <- function(libname, pkgname) {
  with(canprot, {
    # Read amino acid compositions of human proteins and show some information
    human.base <- readRDS(system.file("/extdata/human/human.base.rds", package = "canprot"))
    #packageStartupMessage(paste("human.base:", nrow(human.base), "proteins"))

    human.additional <- readRDS(system.file("/extdata/human/human.additional.rds", package = "canprot"))
    #packageStartupMessage(paste("human.additional:", nrow(human.additional), "proteins"))

    human.extra <- read.csv(system.file("/extdata/human/human.extra.csv", package = "canprot"), as.is = TRUE)
    #packageStartupMessage(paste("human.extra:", nrow(human.extra), "proteins"))

    # Create the data frame with all proteins
    human.aa <- rbind(human.base, human.additional, human.extra)

    # Warn if there are duplicated proteins
    local({
      idup <- duplicated(human.aa$protein)
      if(any(idup)) warning("canprot/zzz.R: duplicated proteins in human.aa: ", paste(human.aa$protein[idup], collapse = " "))
    })
  })
}
