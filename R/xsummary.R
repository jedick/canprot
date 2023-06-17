# canprot/R/xsummary.R
# Format the summary table using xtable
# 20160709 jmd

xsummary <- function(comptab, vars=c("Zc", "nH2O")) {
  # Convert to data frame if needed
  if(!is.data.frame(comptab)) comptab <- do.call(rbind, comptab)
  # Create letter labels
  rownames(comptab) <- c(letters, LETTERS)[1:nrow(comptab)]
  # Return this data frame when the function exits
  comptab.out <- comptab
  # Get the publication key from the dataset name
  publication <- sapply(strsplit(comptab$dataset, "_"), "[", 1)
  # Format the publication key (monospaced font)
  publication <- paste0("<code>", publication, "</code>")
  # Combine the publication and description
  comptab$description <- paste0(publication, " (", comptab$description, ")")
  # Datasets have letter labels
  comptab$dataset <- c(letters, LETTERS)[1:nrow(comptab)]
  # To save column width, change "dataset" to "set"
  colnames(comptab)[1] <- "set"
  # Select the columns to print
  comptab <- comptab[, c(1:4, 7:9, 12:14)]
  # Round values in some columns
  comptab[, c(5, 8)] <- round(comptab[, c(5, 8)], 3)    # mean difference
  comptab[, c(6, 9)] <- signif(comptab[, c(6, 9)], 2)   # CLES
  # Place a marker around high effect size and low p-value
  for(k in vars) {
    # Effect size
    jes <- paste0(k, ".CLES")
    ihigh <- abs(comptab[, jes] - 50) >= 10
    comptab[ihigh, jes] <- paste("**", comptab[ihigh, jes], "**")
    # p-value
    jpv <- paste0(k, ".p.value")
    ilow <- comptab[, jpv] < 0.05
    comptab[, jpv] <- format(comptab[, jpv], digits=1)
    comptab[ilow, jpv] <- paste("**", comptab[ilow, jpv], "**")
    # Bold mean difference if both effect size and p-value are highlighted
    jmd <- paste0(k, ".diff")
    comptab[, jmd] <- format(comptab[, jmd])
    comptab[ihigh & ilow, jmd] <- paste("**", comptab[ihigh & ilow, jmd], "**") 
    # Underline mean difference if only p-value or only effect size is highlighted
    comptab[xor(ilow, ihigh), jmd] <- paste("++", comptab[xor(ilow, ihigh), jmd], "++") 
  }
  # Create xtable
  x <- xtable::xtable(comptab, align=c("c", "l", "l", rep("r", 8)))
  x <- capture.output(xtable::print.xtable(x, type="html",
                                   include.rownames=FALSE,
                                   math.style.exponents=TRUE,
                                   sanitize.text.function=function(x){x}))
  # Bold the indicated effect size, p-values and mean differences
  x <- gsub("> ** ", "> <B>", x, fixed=TRUE)
  x <- gsub(" ** <", "</B> <", x, fixed=TRUE)
  # Underline the indicated mean differences
  x <- gsub("> ++ ", "> <U>", x, fixed=TRUE)
  x <- gsub(" ++ <", "</U> <", x, fixed=TRUE)
  # Add headers that span multiple columns
  span_empty <- "<th colspan=\"4\"></th>"
  span_Zc <- "<th colspan=\"3\"><i>Z</i><sub>C</sub></th>"
  span_nH2O <- "<th colspan=\"3\"><i>n</i><sub>H<sub>2</sub>O</sub></sub></th>"
  span_nO2 <- "<th colspan=\"3\"><i>n</i><sub>O<sub>2</sub></sub></sub></th>"
  span_nAA <- "<th colspan=\"3\"><i>n</i><sub>AA</sub></th>"
  span_var1 <- get(paste0("span_", vars[1]))
  span_var2 <- get(paste0("span_", vars[2]))
  x <- gsub("<table border=1>",
            paste("<table border=1> <tr>", span_empty, span_var1, span_var2, "</tr>"),
            x, fixed=TRUE)
  # More formatting of the headers
  x <- gsub("description", "reference (description)", x, fixed=TRUE)
  x <- gsub("n1", "<i>n</i><sub>1</sub>", x, fixed=TRUE)
  x <- gsub("n2", "<i>n</i><sub>2</sub>", x, fixed=TRUE)
  x <- gsub(paste0(vars[1], ".diff"), "MD", x, fixed=TRUE)
  x <- gsub(paste0(vars[2], ".diff"), "MD", x, fixed=TRUE)
  x <- gsub(paste0(vars[1], ".CLES"), "ES", x, fixed=TRUE)
  x <- gsub(paste0(vars[2], ".CLES"), "ES", x, fixed=TRUE)
  x <- gsub(paste0(vars[1], ".p.value"), "<i>p</i>-value", x, fixed=TRUE)
  x <- gsub(paste0(vars[2], ".p.value"), "<i>p</i>-value", x, fixed=TRUE)
  # Take out extraneous spaces (triggers pre-formatted text - in markdown?)
  x <- gsub("   ", " ", x, fixed=TRUE)
  # Done!
  write.table(x, "", row.names=FALSE, col.names=FALSE, quote=FALSE)
  return(invisible(comptab.out))
}
