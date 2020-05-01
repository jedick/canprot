# canH2O/R/xsummary2.R
# make table for new vignettes
# 20191206

xsummary2 <- function(comptab1, comptab2, comptab3, comptab4) {
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

    PS_TPPG17.down = ct3$PS_TPPG17.mean1,
    PS_TPPG17.up = ct3$PS_TPPG17.mean2,
    PS_TPPG17.diff = ct3$PS_TPPG17.diff,

    PS_LMM16.down = ct3$PS_LMM16.mean1,
    PS_LMM16.up = ct3$PS_LMM16.mean2,
    PS_LMM16.diff = ct3$PS_LMM16.diff,

    nAA.down = ct4$nAA.median1,
    nAA.up = ct4$nAA.median2,
    nAA.diff = ct4$nAA.diff,

    MW.down = ct4$MW.median1,
    MW.up = ct4$MW.median2,
    MW.diff = ct4$MW.diff,

    stringsAsFactors = FALSE
  )

  # prepare table
  x <- out[, c("dataset", "description", "n1", "n2", "ZC.diff", "nH2O_rQEC.diff",
               "nO2_biosynth.diff", "nH2O_biosynth.diff", "nAA.diff", "PS_TPPG17.diff", "PS_LMM16.diff")]
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
  # multiply values of ZC, nO2, nH2O by 1000
  x[, 5:8] <- x[, 5:8] * 1000
  # multiply values of PS by 100
  x[, 10:11] <- x[, 10:11] * 100
  # round values
  x[, 5:11] <- round(x[, 5:11])

  # put markers around negative values
  for(icol in 5:11) {
    ineg <- x[, icol] < 0
    ineg[is.na(ineg)] <- FALSE
    x[ineg, icol] <- paste("**", x[ineg, icol], "**")
  }

  # remove PS columns if all values are NA 20200101
  if(all(is.na(x[, 10:11]))) x <- x[, -(10:11)]

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
  span_empty1 <- "<td align=\"center\" colspan=\"1\"></td>"
  span_TPPG17 <- "<td align=\"center\" colspan=\"1\"><B>TPPG17</B></td>"
  span_LMM16 <- "<td align=\"center\" colspan=\"1\"><B>LMM16</B></td>"
  border <- paste("<table border=1> <tr>", span_empty5, span_rQEC, span_biosynth, span_empty1, span_TPPG17, span_LMM16, "</tr>")
  if(!any(grepl("PS_", x))) border <- paste("<table border=1> <tr>", span_empty5, span_rQEC, span_biosynth, span_empty1, "</tr>")
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
  x <- gsub("PS_TPPG17.diff", "&Delta;PS", x, fixed=TRUE)
  x <- gsub("PS_LMM16.diff", "&Delta;PS", x, fixed=TRUE)

  # done!
  write.table(x, "", row.names=FALSE, col.names=FALSE, quote=FALSE)
  return(invisible(out))
}
