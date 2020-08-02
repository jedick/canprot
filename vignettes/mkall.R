# canprot/vignettes/mkall.R
# Process vignette and move html file to ../inst/doc
# 20200802

files <- dir(pattern = "Rmd")
print(system.time(
  for(f in files) {
    sep <- paste(rep("=", nchar(f) + 4), collapse = "")
    message()
    print(sep, quote = FALSE)
    print(paste("|", f, "|"), quote = FALSE)
    print(sep, quote = FALSE)
    rmarkdown::render(f)
    from <- gsub("Rmd", "html", f)
    to <- file.path("../inst/doc", from)
    file.rename(from, to)
  }
))
