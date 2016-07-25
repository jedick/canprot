# canprot/R/CLES.R
# calculate common language effect size
# 20160708 jmd

CLES <- function(x, y) {
  # generate all combinations of elements in x and y
  combinations <- expand.grid(x=x, y=y)
  # calculate the differences
  differences <- combinations$y - combinations$x
  # find the proportion of positive differences
  proportion <- sum(differences > 0) / length(differences)
  return(proportion)
}
