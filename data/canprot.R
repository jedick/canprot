# canprot/data/canprot.R
# create the data objects for amino acid compositions of human proteins
# 20160705 jmd

# make sure the CHNOSZ environment is attached
if(!"CHNOSZ" %in% search()) {
  data(thermo, package="CHNOSZ")
}

# create the canprot environment if it does not exist
if(!"canprot" %in% search()) {
  attach(NULL, name="canprot")
  message("data(canprot): attached environment \"canprot\"")
  with(as.environment("canprot"), {
    # read amino acid compositions of human proteins and show some information
    #human_base <- read.csv("human_base.csv.xz", as.is=TRUE)
    load("human_base.Rdata")
    message(paste("human_base:", nrow(human_base), "proteins"))
    #human_additional <- read.csv("human_additional.csv.xz", as.is=TRUE)
    load("human_additional.Rdata")
    message(paste("human_additional:", nrow(human_additional), "proteins"))
    human_extra <- read.csv("human_extra.csv", as.is=TRUE)
    message(paste("human_extra:", nrow(human_extra), "proteins"))
    # create the object
    human_aa <- rbind(human_base, human_additional, human_extra)
    # warn if there are duplicated proteins
    local({
      idup <- duplicated(human_aa$protein)
      if(any(idup)) warning("data(human): duplicated proteins: ", paste(human_aa$protein[idup], collapse=" "))
    })
    # updates for UniProt IDs
    uniprot_updates <- read.csv("uniprot_updates.csv", as.is=TRUE)
  })
}
