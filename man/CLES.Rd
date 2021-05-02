\encoding{UTF-8}
\name{CLES}
\alias{CLES}
\title{Common Language Effect Size}
\description{
  Calculate the common language effect size.
}

\usage{
  CLES(x, y, distribution = "normal")
}

\arguments{
  \item{x}{numeric, data}
  \item{y}{numeric, data}
  \item{distribution}{\samp{normal} to use probabilities calculated for a normal distribution, or NA for empirical probabilities}
}

\details{
The common language statistic is defined for continuous data as \dQuote{the probability that a score sampled at random from one distribution will be greater than a score sampled from some other distribution} (McGraw and Wong, 1992).

Given the default value of \code{distribution} (\samp{normal}), this function uses \code{\link{pnorm}} to calculate the probability that a random sample from the unit normal distribution is greater than the Z score (i.e. (the mean of \samp{y} minus the mean of \samp{x}) / square root of (variance of \samp{x} plus variance of \samp{y})).

If \code{distribution} is NA, this function calculates the empirical probability that the difference is positive, that is, the fraction of all possible pairings between elements of \code{x} and \code{y} where the difference (\samp{y} value - \samp{x} value) is positive.
It may not be possible to calculate the empirical probability for very large samples because of memory limits.

The examples use \emph{simulated data for normal distributions}, given the sample size, mean, and standard deviation of datasets cited by McGraw and Wong, 1992.
Therefore, the empirical probability in the examples approaches the normal curve probability.
However, the empirical probability for \emph{nonnormal} distributions is distinct from the normal curve probability, as discussed on p. 364-365 of McGraw and Wong, 1992.
}

\examples{
# Example 1: Height differences between males and females
# a) Use statistics quoted by McGraw and Wong, 1992 from NCHS, 1987
# for heights in inches of 18-24 year-old males and females
# Table 14: number, mean height, and standard deviation of height of females
n1 <- 1066
M1 <- 64.3
SD1 <- 2.8
# Table 13: number, mean height, and standard deviation of height of males
n2 <- 988
M2 <- 69.7
SD2 <- 2.6
# b) Simulate data from a normal distribution with exact mean and SD
# use rnorm2 function from Ben Bolker's answer to
# https://stackoverflow.com/questions/18919091/generate-random-numbers-with-fixed-mean-and-sd
rnorm2 <- function(n, mean, sd) { mean + sd * scale(rnorm(n)) }
set.seed(1234)
height_female <- rnorm2(n1, M1, SD1)
height_male <- rnorm2(n2, M2, SD2)
# c) Calculate the CLES using the normal distribution and empirical probability
CLES_normal <- CLES(height_female, height_male)
CLES_empirical <- CLES(height_female, height_male, distribution = NA)
# d) Test numerical equivalence of the results
# The CLES is approximately 0.92 (McGraw and Wong, 1992)
# (note: becasue we used rnorm2, this doesn't depend on the seed)
stopifnot(all.equal(CLES_normal, 0.92, tol = 0.01))
# With this seed, the difference between the normal curve probability
# and empirical probability is less than 1%
stopifnot(all.equal(CLES_normal, CLES_empirical, tol = 0.01))

# Example 1.5: Use multiple simulated datasets to show approach
# of empirical probability to normal curve probability
CLES_empirical_n <- sapply(1:100, function(x) {
  height_female <- rnorm2(n1, M1, SD1)
  height_male <- rnorm2(n2, M2, SD2)
  CLES(height_female, height_male, distribution = NA)
})
CLES_empirical <- mean(CLES_empirical_n)
# now we're even closer to the normal curve probability
stopifnot(all.equal(CLES_normal, CLES_empirical, tol = 0.0001))

# Example 2: Multiple datasets in Table 2 of McGraw and Wong, 1992
# Sample statistics for females
n1 <- c(638, 672, 3139, 420740, 19274, 104263, 207, 394, 1066, 982, 108, 108)
M1 <- c(103, 15, 103, 18.9, 30, 16.1, 6.9, 13.3, 64.3, 134, 45, 94)
SD1 <- sqrt(c(908, 74, 219, 27, 110, 59, 15, 164, 6.8, 688, 310, 1971))
# Sample statistics for males
n2 <- c(354, 359, 3028, 356704, 21768, 133882, 199, 469, 988, 988, 443, 443)
M2 <- c(112, 23, 100, 17.9, 33, 18.6, 9.3, 21.8, 69.7, 163, 86, 212)
SD2 <- sqrt(c(1096, 96, 202, 29, 110, 61, 15, 133, 7.8, 784, 818, 5852))
# A function to calculate the effect size using simulated data
CLESfun <- function(n1, M1, SD1, n2, M2, SD2, distribution) {
  rnorm2 <- function(n, mean, sd) { mean + sd * scale(rnorm(n)) }
  set.seed(1234)
  x <- rnorm2(n1, M1, SD1)
  y <- rnorm2(n2, M2, SD2)
  CLES(x, y, distribution)
}
# Calculate 100 * CL for the normal curve probabilities
CLnorm <- sapply(1:12, function(i) {
  CL <- CLESfun(n1[i], M1[i], SD1[i], n2[i], M2[i], SD2[i], "normal")
  round(100 * CL)
})
# Calculate 100 * CL for empirical probabilities
CLemp <- sapply(1:12, function(i) {
  # skip very large samples: not enough memory
  if(n1[i] > 5000 | n2[i] > 5000) NA else {
    CL <- CLESfun(n1[i], M1[i], SD1[i], n2[i], M2[i], SD2[i], NA)
    round(100 * CL)
  }
})
# The difference between the empirical and normal curve
# probabilities is not more than 1 percent
stopifnot(max(abs(CLemp - CLnorm), na.rm = TRUE) <= 1)
# TODO: Why are some of the calculated values different from
# Table 2 of McGraw and Wong, 1992?
CLref <- c(54, 74, 44, 45, 56, 63, 67, 65, 92, 78, 89, 91)
# Differences range from -4 to 4
range(CLnorm - CLref) 
#stopifnot(max(abs(CLnorm - CLref)) == 0)
}

\references{
  McGraw, Kenneth O. and Wong, S. P. (1992) A common language effect size statistic. \emph{Psychological Bulletin} \bold{11}, 361--365. \doi{10.1037/0033-2909.111.2.361}

  National Center for Health Statistics (1987) \emph{Anthropometric Reference Data and Prevalence of Overweight: United States, 1976-1980.} Data from the National Health Survey, Series 11, No. 238. DHHS Publication (PHS) No. 87-1688. U.S. Government Printing Office, Washington, DC. \url{https://www.cdc.gov/nchs/data/series/sr_11/sr11_238.pdf}

}

\concept{Statistics}
