# canH2O/R/xsummary3.R
# make table for new vignettes (with GRAVY and pI)
# 20200418

xsummary3 <- function(comptab1, comptab2, comptab3) {
  ct1 <- do.call(rbind, comptab1)
  ct2 <- do.call(rbind, comptab2)
  ct3 <- do.call(rbind, comptab3)
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

    pI.down = ct2$pI.median1,
    pI.up = ct2$pI.median2,
    pI.diff = ct2$pI.diff,

    GRAVY.down = ct2$GRAVY.median1,
    GRAVY.up = ct2$GRAVY.median2,
    GRAVY.diff = ct2$GRAVY.diff,

    nAA.down = ct3$nAA.median1,
    nAA.up = ct3$nAA.median2,
    nAA.diff = ct3$nAA.diff,

    MW.down = ct3$MW.median1,
    MW.up = ct3$MW.median2,
    MW.diff = ct3$MW.diff,

    stringsAsFactors = FALSE
  )

  # prepare table
  x <- out[, c("dataset", "description", "n1", "n2", "ZC.diff", "nH2O_rQEC.diff",
               "pI.diff", "GRAVY.diff", "nAA.diff", "MW.diff")]
  # get the publication key from the dataset name
  publication <- sapply(strsplit(x$dataset, "_"), "[", 1)
  # format the publication key (monospaced font)
  publication <- paste0("<code>", publication, "</code>")
  # italicize species names (first words in description, surrounded by underscore)
  x$description <- gsub("^_", "<i>", x$description)
  x$description <- gsub("_", "</i>", x$description)
  # combine the publication and description
  x$description <- paste0(publication, " (", x$description, ")")
  # datasets have letter labels
  x$dataset <- c(letters, LETTERS)[1:nrow(x)]
  # to save column width, change "dataset" to "set"
  colnames(x)[1] <- "set"
  # multiply values of ZC, nH2O and GRAVY by 1000
  x[, c(5:6, 8)] <- x[, c(5:6, 8)] * 1000
  # multiply values of pI and MW by 100
  x[, c(7, 10)] <- x[, c(7, 10)] * 100
  # round values
  x[, 5:10] <- round(x[, 5:10])

  # put markers around negative values
  for(icol in 5:10) {
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

#  # add headers that span multiple columns
#  span_empty5 <- "<td align=\"center\" colspan=\"5\"></td>"
#  span_rQEC <- "<td align=\"center\" colspan=\"1\"><B>rQEC</B></td>"
#  span_biosynth <- "<td align=\"center\" colspan=\"2\"><B>biosynthetic</B></td>"
#  span_empty3 <- "<td align=\"center\" colspan=\"3\"></td>"
#  border <- paste("<table border=1> <tr>", span_empty5, span_rQEC, span_biosynth, span_empty3, "</tr>")
#  x <- gsub("<table border=1>", border, x, fixed=TRUE)

  # more formatting of the headers
  x <- gsub("description", "reference (description)", x, fixed=TRUE)
  x <- gsub("n1", "<i>n</i><sub>down</sub>", x, fixed=TRUE)
  x <- gsub("n2", "<i>n</i><sub>up</sub>", x, fixed=TRUE)
  x <- gsub("ZC.diff", "&Delta;<i>Z</i><sub>C</sub>", x, fixed=TRUE)
  x <- gsub("nH2O_rQEC.diff", "&Delta;<i>n</i><sub>H<sub>2</sub>O</sub>", x, fixed=TRUE)
  x <- gsub("pI.diff", "&Delta;pI", x, fixed=TRUE)
  x <- gsub("GRAVY.diff", "&Delta;GRAVY", x, fixed=TRUE)
  x <- gsub("nAA.diff", "&Delta;<i>n</i><sub>AA</sub>", x, fixed=TRUE)
  x <- gsub("nAA.diff", "&Delta;MW", x, fixed=TRUE)

  # done!
  write.table(x, "", row.names=FALSE, col.names=FALSE, quote=FALSE)
  return(invisible(out))
}
