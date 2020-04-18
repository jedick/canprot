# canH2O/R/xsummary3.R
# make table for new vignettes (with GRAVY and pI)
# 20200418

xsummary3 <- function(comptab1, comptab2, comptab3, comptab4) {
  ct1 <- do.call(rbind, comptab1)
  ct2 <- do.call(rbind, comptab2)
  ct3 <- do.call(rbind, comptab3)
  ct4 <- do.call(rbind, comptab4)
  # get all data
  # include medians or means for up and down groups for making summary .csv files 20200125
  out <- data.frame(
    dataset = ct1$dataset,
    description = ct2$description,
    n1 = ct1$n1,
    n2 = ct1$n2,

    ZC.down = ct1$ZC.median1,
    ZC.up = ct1$ZC.median2,
    ZC.diff = ct1$ZC.diff,

    nH2O_rQEC.down = ct1$nH2O.median1,
    nH2O_rQEC.up = ct1$nH2O.median2,
    nH2O_rQEC.diff = ct1$nH2O.diff,

    nO2_biosynth.down = ct2$nO2.median1,
    nO2_biosynth.up = ct2$nO2.median2,
    nO2_biosynth.diff = ct2$nO2.diff,

    nH2O_biosynth.down = ct2$nH2O.median1,
    nH2O_biosynth.up = ct2$nH2O.median2,
    nH2O_biosynth.diff = ct2$nH2O.diff,

    nAA.down = ct3$nAA.median1,
    nAA.up = ct3$nAA.median2,
    nAA.diff = ct3$nAA.diff,

    pI.down = ct4$pI.median1,
    pI.up = ct4$pI.median2,
    pI.diff = ct4$pI.diff,

    GRAVY.down = ct4$GRAVY.median1,
    GRAVY.up = ct4$GRAVY.median2,
    GRAVY.diff = ct4$GRAVY.diff,

    stringsAsFactors = FALSE
  )

  # prepare table
  x <- out[, c("dataset", "description", "n1", "n2", "ZC.diff", "nH2O_rQEC.diff",
               "nO2_biosynth.diff", "nH2O_biosynth.diff", "nAA.diff", "pI.diff", "GRAVY.diff")]
  # get the publication key from the dataset name
  publication <- sapply(strsplit(x$dataset, "_"), "[", 1)
  # format the publication key (monospaced font)
  publication <- paste0("<code>", publication, "</code>")
  # combine the publication and description
  x$description <- paste0(publication, " (", x$description, ")")
  # datasets have letter labels
  x$dataset <- c(letters, LETTERS)[1:nrow(x)]
  # to save column width, change "dataset" to "set"
  colnames(x)[1] <- "set"
  # multiply values of ZC, nO2, nH2O by 1000
  x[, 5:8] <- x[, 5:8] * 1000
  # multiply values of pI by 100
  x[, 10] <- x[, 10] * 100
  # multiply values of GRAVY by 1000
  x[, 11] <- x[, 11] * 1000
  # round values
  x[, 5:11] <- round(x[, 5:11])

  # put markers around negative values
  for(icol in 5:11) {
    ineg <- x[, icol] < 0
    ineg[is.na(ineg)] <- FALSE
    x[ineg, icol] <- paste("**", x[ineg, icol], "**")
  }

  # create xtable
  x <- xtable::xtable(x, align=c("c", "l", "l", rep("r", ncol(x) - 2)))
  x <- capture.output(xtable::print.xtable(x, type="html",
                                   include.rownames=FALSE,
                                   math.style.exponents=TRUE,
                                   sanitize.text.function=function(x){x}))

  # make the marked negative values bold
  x <- gsub("> ** ", "> <B>", x, fixed=TRUE)
  x <- gsub(" ** <", "</B> <", x, fixed=TRUE)
  # change "NaN" to "NA"
  x <- gsub("NaN", "NA", x, fixed=TRUE)

  # add headers that span multiple columns
  span_empty5 <- "<td align=\"center\" colspan=\"5\"></td>"
  span_rQEC <- "<td align=\"center\" colspan=\"1\"><B>rQEC</B></td>"
  span_biosynth <- "<td align=\"center\" colspan=\"2\"><B>biosynthetic</B></td>"
  span_empty3 <- "<td align=\"center\" colspan=\"3\"></td>"
  border <- paste("<table border=1> <tr>", span_empty5, span_rQEC, span_biosynth, span_empty3, "</tr>")
  x <- gsub("<table border=1>", border, x, fixed=TRUE)

  # more formatting of the headers
  x <- gsub("description", "reference (description)", x, fixed=TRUE)
  x <- gsub("n1", "<i>n</i><sub>down</sub>", x, fixed=TRUE)
  x <- gsub("n2", "<i>n</i><sub>up</sub>", x, fixed=TRUE)
  x <- gsub("ZC.diff", "&Delta;<i>Z</i><sub>C</sub>", x, fixed=TRUE)
  x <- gsub("nH2O_rQEC.diff", "&Delta;<i>n</i><sub>H<sub>2</sub>O</sub>", x, fixed=TRUE)
  x <- gsub("nO2_biosynth.diff", "&Delta;<i>n</i><sub>O<sub>2</sub></sub>", x, fixed=TRUE)
  x <- gsub("nH2O_biosynth.diff", "&Delta;<i>n</i><sub>H<sub>2</sub>O</sub>", x, fixed=TRUE)
  x <- gsub("nAA.diff", "&Delta;<i>n</i><sub>AA</sub>", x, fixed=TRUE)
  x <- gsub("pI.diff", "&Delta;pI", x, fixed=TRUE)
  x <- gsub("GRAVY.diff", "&Delta;GRAVY", x, fixed=TRUE)

  # done!
  write.table(x, "", row.names=FALSE, col.names=FALSE, quote=FALSE)
  return(invisible(out))
}
