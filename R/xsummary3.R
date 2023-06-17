# canH2O/R/xsummary3.R
# Make table for new vignettes (with GRAVY and pI)
# 20200418

xsummary3 <- function(comptab1, comptab2, comptab3) {
  ct1 <- do.call(rbind, comptab1)
  ct2 <- do.call(rbind, comptab2)
  ct3 <- do.call(rbind, comptab3)
  # Get all data
  # Include medians or means for up and down groups for making summary .csv files 20200125
  out <- data.frame(
    dataset = ct1$dataset,
    description = ct2$description,
    n1 = ct1$n1,
    n2 = ct1$n2,

    Zc.down = ct1$Zc.median1,
    Zc.up = ct1$Zc.median2,
    Zc.diff = ct1$Zc.diff,

    nH2O.down = ct1$nH2O.median1,
    nH2O.up = ct1$nH2O.median2,
    nH2O.diff = ct1$nH2O.diff,

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

  # Prepare table
  x <- out[, c("dataset", "description", "n1", "n2", "Zc.diff", "nH2O.diff",
               "pI.diff", "GRAVY.diff", "nAA.diff", "MW.diff")]
  # Get the publication key from the dataset name
  publication <- sapply(strsplit(x$dataset, "_"), "[", 1)
  # Format the publication key (monospaced font)
  publication <- paste0("<code>", publication, "</code>")
  # Italicize species names (first words in description, surrounded by underscore)
  x$description <- gsub("^_", "<i>", x$description)
  x$description <- gsub("_", "</i>", x$description)
  # Combine the publication and description
  x$description <- paste0(publication, " (", x$description, ")")
  # Datasets have letter labels
  x$dataset <- c(letters, LETTERS)[1:nrow(x)]
  # To save column width, change "dataset" to "set"
  colnames(x)[1] <- "set"
  # Multiply values of Zc, nH2O and GRAVY by 1000
  x[, c(5:6, 8)] <- x[, c(5:6, 8)] * 1000
  # Multiply values of pI and MW by 100
  x[, c(7, 10)] <- x[, c(7, 10)] * 100
  # Round values
  x[, 5:10] <- round(x[, 5:10])

  # Put markers around negative values
  for(icol in 5:10) {
    ineg <- x[, icol] < 0
    ineg[is.na(ineg)] <- FALSE
    x[ineg, icol] <- paste("**", x[ineg, icol], "**")
  }

  # Create xtable
  x <- xtable::xtable(x, align=c("c", "l", "l", rep("r", ncol(x) - 2)))
  x <- capture.output(xtable::print.xtable(x, type="html",
                                   include.rownames=FALSE,
                                   math.style.exponents=TRUE,
                                   sanitize.text.function=function(x){x}))

  # Make the marked negative values bold
  x <- gsub("> ** ", "> <B>", x, fixed=TRUE)
  x <- gsub(" ** <", "</B> <", x, fixed=TRUE)
  # Change "NaN" to "NA"
  x <- gsub("NaN", "NA", x, fixed=TRUE)

  # More formatting of the headers
  x <- gsub("description", "reference (description)", x, fixed=TRUE)
  x <- gsub("n1", "<i>n</i><sub>down</sub>", x, fixed=TRUE)
  x <- gsub("n2", "<i>n</i><sub>up</sub>", x, fixed=TRUE)
  x <- gsub("Zc.diff", "&Delta;<i>Z</i><sub>C</sub>", x, fixed=TRUE)
  x <- gsub("nH2O.diff", "&Delta;<i>n</i><sub>H<sub>2</sub>O</sub>", x, fixed=TRUE)
  x <- gsub("pI.diff", "&Delta;pI", x, fixed=TRUE)
  x <- gsub("GRAVY.diff", "&Delta;GRAVY", x, fixed=TRUE)
  x <- gsub("nAA.diff", "&Delta;<i>n</i><sub>AA</sub>", x, fixed=TRUE)
  x <- gsub("MW.diff", "&Delta;MW", x, fixed=TRUE)

  # Done!
  write.table(x, "", row.names=FALSE, col.names=FALSE, quote=FALSE)
  return(invisible(out))
}
