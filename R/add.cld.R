# Add compact letter display to boxplot 20220531
# Exported in JMDplots 20220609
# Moved from JMDplots to canprot 20240301
add.cld <- function(datlist, bp, dx = NULL, dy = NULL) {
  # Add names if missing
  if(is.null(names(datlist))) names(datlist) <- 1:length(datlist)
  # Remove empty elements
  datlist_1 <- datlist[sapply(datlist, length) != 0]
  # Turn the list into a data frame with "group" column taken from names of the list elements
  Zcdat <- do.call(rbind, sapply(1:length(datlist_1), function(i) data.frame(group = names(datlist_1)[i], Zc = datlist_1[[i]]), simplify = FALSE))
  # One-way ANOVA and Tukey's Honest Significant Differences
  # Adapted from https://statdoe.com/one-way-anova-and-box-plot-in-r/
  anova <- aov(Zc ~ group, data = Zcdat)
  tukey <- TukeyHSD(anova)
  # Compact letter display
  letters <- multcompLetters4(anova, tukey, reversed = TRUE)$group$Letters
  # Get into same order as data
  letters <- letters[match(names(datlist), names(letters))]
  # Work out dy and dx if NULL
  ngroup <- length(bp$n)
  if(is.null(dx)) dx <- diff(par("usr")[1:2]) / ngroup / 3
  if(is.null(dy)) dy <- diff(par("usr")[3:4]) / 30
  # Add to plot
  text((1:ngroup) + dx, bp$stats[4, ] + dy, letters)
  # Return dx, dy, and letters
  invisible(list(dx = dx, dy = dy, letters = letters))
}
