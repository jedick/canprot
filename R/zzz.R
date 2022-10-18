# canprot/R/zzz.R
# Create the data objects for amino acid compositions of human proteins
# 20160705 jmd (data/canprot.R)
# 20190225 moved to R/zzz.R
# 20200509 change environment name from canprot to human

# The 'human' environment is made here in open code 20190214
# https://stackoverflow.com/questions/41954302/where-to-create-package-environment-variables
human <- new.env()

# Initialize the 'human' environment
.onAttach <- function(libname, pkgname) {
  with(human, {
    # Read amino acid compositions of human proteins and show some information
    human_base <- readRDS(system.file("/extdata/protein/human_base.rds", package = "canprot"))
    packageStartupMessage(paste("human_base:", nrow(human_base), "proteins"))

    human_additional <- readRDS(system.file("/extdata/protein/human_additional.rds", package = "canprot"))
    packageStartupMessage(paste("human_additional:", nrow(human_additional), "proteins"))

    human_extra <- read.csv(system.file("/extdata/protein/human_extra.csv", package = "canprot"), as.is = TRUE)
    packageStartupMessage(paste("human_extra:", nrow(human_extra), "proteins"))

    # create the data frame with all proteins
    human_aa <- rbind(human_base, human_additional, human_extra)

    # Warn if there are duplicated proteins
    local({
      idup <- duplicated(human_aa$protein)
      if(any(idup)) warning("data(human): duplicated proteins: ", paste(human_aa$protein[idup], collapse = " "))
    })

    # Load updates for UniProt IDs
    uniprot_updates <- read.csv(system.file("/extdata/protein/uniprot_updates.csv", package = "canprot"), as.is = TRUE)
  })
}

# Set 'basis' option 20200817
# Adapted from R/src/library/grDevices/zzz.R
.onLoad <- function(libname, pkgname) {
  op.canprot <- list(basis = "QEC")
  toset <- !(names(op.canprot) %in% names(.Options))
  if(any(toset)) options(op.canprot[toset])
}
