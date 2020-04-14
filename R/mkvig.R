# canprot/mkvig.R
# compile and view vignettes from command line
# 20200414 jmd

mkvig <- function(vig = "colorectal") {
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
