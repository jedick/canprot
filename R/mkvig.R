# canprot/mkvig.R
# compile and view vignettes from command line
# 20200414 jmd

mkvig <- function(vig = NULL) {
  vig.allowed <- gsub("pdat_", "", grep("pdat_", ls("package:canprot"), value = TRUE))
  isnull <- is.null(vig)
  toomany <- length(vig) > 1
  notallowed <- !any(vig %in% vig.allowed)
  if(isnull | toomany | notallowed) stop("'vig' should be one of: ", paste(vig.allowed, collapse = ", "))
  # location of the vignette directory
  vigdir <- system.file("vignettes", package = "canprot")
  # names of the vignette source and html output files
  vigfile <- file.path(vigdir, paste0(vig, ".Rmd"))
  htmlfile <- tempfile(pattern = paste0(vig, "_"), fileext = ".html")
  # compile the vignette and open it in the browser
  rmarkdown::render(vigfile, output_file = htmlfile, knit_root_dir = vigdir)
  browseURL(htmlfile)
  # return the path of html file
  htmlfile
}

# mkoldvig added 20200829
mkoldvig <- function(vig = NULL) {
  vig.allowed <- c("colorectal", "pancreatic", "hypoxia", "hyperosmotic")
  isnull <- is.null(vig)
  toomany <- length(vig) > 1
  notallowed <- !any(vig %in% vig.allowed)
  if(isnull | toomany | notallowed) stop("'vig' should be one of: ", paste(vig.allowed, collapse = ", "))
  # location of the vignette directory
  vigdir <- system.file("oldvignettes", package = "canprot")
  # names of the vignette source and html output files
  vigfile <- file.path(vigdir, paste0(vig, ".Rmd"))
  htmlfile <- tempfile(pattern = paste0(vig, "_"), fileext = ".html")
  # compile the vignette and open it in the browser
  rmarkdown::render(vigfile, output_file = htmlfile, knit_root_dir = vigdir)
  browseURL(htmlfile)
  # return the path of html file
  htmlfile
}
