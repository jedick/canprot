# canprot/R/internal.R
# internal (non-exported) objects
# 20160706 jmd

# colors for points on scatterplots
cpcol <- list(
  # colorspace::diverge_hcl(2, c=100, l=c(50, 90), power=1)[1],
  blue = "#4A6FE3",
  # colorspace::diverge_hcl(2, c=100, l=c(50, 90), power=1)[2]
  red = "#D33F6A"
)

# "oldstyle" labels including overbar
cplabbar <- cplab
cplabbar$nH2O <- expression(bar(italic(n))[H[2]*O])
cplabbar$DnH2O <- expression(Delta*bar(italic(n))[H[2]*O])
