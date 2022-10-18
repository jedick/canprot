# canprot/mkvig.R
# Compile and view vignettes from command line
# 20200414 jmd

mkvig <- function(vig = NULL) {
  vig.allowed <- gsub("pdat_", "", grep("pdat_", ls("package:canprot"), value = TRUE))
  isnull <- is.null(vig)
  toomany <- length(vig) > 1
  notallowed <- !any(vig %in% vig.allowed)
  if(isnull | toomany | notallowed) stop("'vig' should be one of: ", paste(vig.allowed, collapse = ", "))
  # Location of the vignette directory
  vigdir <- system.file("vignettes", package = "canprot")
  # Names of the vignette source and html output files
  vigfile <- file.path(vigdir, paste0(vig, ".Rmd"))
  htmlfile <- tempfile(pattern = paste0(vig, "_"), fileext = ".html")
  # Compile the vignette and open it in the browser
  rmarkdown::render(vigfile, output_file = htmlfile, knit_root_dir = vigdir)
  # Set 'browser' option so vignettes open under Rstudio 20201017
  # https://stackoverflow.com/questions/62536479/the-command-exams2html-does-not-generate-html-page-when-it-is-run-from-rstudio
  if(.Platform$OS.type == "windows") oldopt <- options(browser = NULL)
  browseURL(htmlfile)
  if(.Platform$OS.type == "windows") options(oldopt)
  # Return the path of html file
  htmlfile
}
