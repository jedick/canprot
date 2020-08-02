# canH2O/vignettes/mkall.R
# use this to process all vignettes in order to make
# the csv files that are installed with the package

files <- dir(pattern = "Rmd")
print(system.time(
  for(f in files) {
    sep <- paste(rep("=", nchar(f) + 4), collapse = "")
    message()
    print(sep, quote = FALSE)
    print(paste("|", f, "|"), quote = FALSE)
    print(sep, quote = FALSE)
    rmarkdown::render(f)
  }
))
