# canprot/R/CLES.R
# calculate common language effect size
# 20160708 jmd first version, empirical probability only
# 20200528 add 'distribution' argument and normal curve probability

CLES <- function(x, y, distribution = "normal") {
  if(identical(distribution, "normal")) {
    ## Normal curve probability
    # Calculate mean and SD of samples
    M1 <- mean(x)
    M2 <- mean(y)
    SD1 <- sd(x)
    SD2 <- sd(y)
    # Calculate mean and SD of differences
    # (e.g. center plot of McGraw and Wong, 1992, Figure 1)
    M <- M2 - M1
    SD <- sqrt(SD1^2 + SD2^2)
    # Calculate Z score
    Z <- M / SD
    # Calculate probability
    pnorm(Z)
  } else if(identical(distribution, NA)) {
    ## Empirical probability
    # generate all combinations of elements in x and y
    combinations <- expand.grid(x = x, y = y)
    # calculate the differences
    differences <- combinations$y - combinations$x
    # find the proportion of positive differences
    sum(differences > 0) / length(differences)
  } else stop("distribution should be 'normal' or NA")
}
