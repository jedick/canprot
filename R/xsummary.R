# canprot/R/xsummary.R
# format the summary table using xtable
# 20160709 jmd

xsummary <- function(comptab) {
  # convert to data frame if needed
  if(!is.data.frame(comptab)) comptab <- do.call(rbind, comptab)
  # create letter labels
  rownames(comptab) <- c(letters, LETTERS)[1:nrow(comptab)]
  # return this data frame when the function exits
  comptab.out <- comptab
  # get the publication key from the dataset name
  publication <- sapply(strsplit(comptab$dataset, "_"), "[", 1)
  # format the publication key (monospaced font)
  publication <- paste0("<code>", publication, "</code>")
  # combine the publication and description
  comptab$description <- paste0(publication, " (", comptab$description, ")")
  # datasets have letter labels
  comptab$dataset <- c(letters, LETTERS)[1:nrow(comptab)]
  # to save column width, change "dataset" to "set"
  colnames(comptab)[1] <- "set"
  # select the columns to print
  comptab <- comptab[, c(1:4, 7:9, 12:14)]
  # round values in some columns
  comptab[, c(5, 8)] <- round(comptab[, c(5, 8)], 3)    # mean difference
  comptab[, c(6, 9)] <- signif(comptab[, c(6, 9)], 2)   # CLES
  # place a marker around high effect size and low p-value
  for(k in c("ZC", "nH2O")) {
    # effect size
    jes <- paste0(k, ".CLES")
    ihigh <- abs(comptab[, jes] - 50) >= 10
    comptab[ihigh, jes] <- paste("**", comptab[ihigh, jes], "**")
    # p-value
    jpv <- paste0(k, ".p.value")
    ilow <- comptab[, jpv] < 0.05
    comptab[, jpv] <- format(comptab[, jpv], digits=1)
    comptab[ilow, jpv] <- paste("**", comptab[ilow, jpv], "**")
    # bold mean difference if both effect size and p-value are highlighted
    jmd <- paste0(k, ".diff")
    comptab[, jmd] <- format(comptab[, jmd])
    comptab[ihigh & ilow, jmd] <- paste("**", comptab[ihigh & ilow, jmd], "**") 
    # underline mean difference if only p-value or only effect size is highlighted
    comptab[xor(ilow, ihigh), jmd] <- paste("++", comptab[xor(ilow, ihigh), jmd], "++") 
  }
  # create xtable
  x <- xtable::xtable(comptab, align=c("c", "l", "l", rep("r", 8)))
  x <- capture.output(xtable::print.xtable(x, type="html",
                                   include.rownames=FALSE,
                                   math.style.exponents=TRUE,
                                   sanitize.text.function=function(x){x}))
  # bold the indicated effect size, p-values and mean differences
  x <- gsub("> ** ", "> <B>", x, fixed=TRUE)
  x <- gsub(" ** <", "</B> <", x, fixed=TRUE)
  # underling the indicated mean differences
  x <- gsub("> ++ ", "> <U>", x, fixed=TRUE)
  x <- gsub(" ++ <", "</U> <", x, fixed=TRUE)
  # add headers that span multiple columns
  span_empty <- "<th colspan=\"4\"></th>"
  span_ZC <- "<th colspan=\"3\"><i>Z<i><sub>C</sub></th>"
  span_nH2O <- "<th colspan=\"3\"><i>n<i><sub>H<sub>2</sub>O</sub></sub></th>"
  x <- gsub("<table border=1>",
            paste("<table border=1> <tr>", span_empty, span_ZC, span_nH2O, "</tr>"),
            x, fixed=TRUE)
  # more formatting of the headers
  x <- gsub("description", "reference (description)", x, fixed=TRUE)
  x <- gsub("n1", "<i>n</i><sub>1</sub>", x, fixed=TRUE)
  x <- gsub("n2", "<i>n</i><sub>2</sub>", x, fixed=TRUE)
  x <- gsub("ZC.diff", "MD", x, fixed=TRUE)
  x <- gsub("nH2O.diff", "MD", x, fixed=TRUE)
  x <- gsub("ZC.CLES", "ES", x, fixed=TRUE)
  x <- gsub("nH2O.CLES", "ES", x, fixed=TRUE)
  x <- gsub("ZC.p.value", "<i>p</i>-value", x, fixed=TRUE)
  x <- gsub("nH2O.p.value", "<i>p</i>-value", x, fixed=TRUE)
  # done!
  write.table(x, "", row.names=FALSE, col.names=FALSE, quote=FALSE)
  return(invisible(comptab.out))
}
