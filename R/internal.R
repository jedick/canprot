# canprot/R/internal.R
# internal (non-exported) objects
# 20160706 jmd

# colors for points on scatterplots
cpcol <- list(
  blue = colorspace::diverge_hcl(2, c=100, l=c(50, 90), power=1)[1],
  red = colorspace::diverge_hcl(2, c=100, l=c(50, 90), power=1)[2]
)

# text for figure labels
cplab <- list(
  nH2O = expression(bar(italic(n))[H[2]*O]),
  DnH2O = expression(Delta*bar(italic(n))[H[2]*O]),
  ZC = expression(italic(Z)[C]),
  DZC = expression(Delta*italic(Z)[C]),
  logfO2 = expression(log*italic("f")[O[2]*group("(", italic("g"), ")")]),
  logaH2O = expression(log*italic("a")[H[2]*O*group("(", italic("liq"), ")")]),
  nC = expression(bar(italic(n))[C]),
  nN = expression(bar(italic(n))[N]),
  nS = expression(bar(italic(n))[S]),
  DnC = expression(Delta*bar(italic(n))[C]),
  DnN = expression(Delta*bar(italic(n))[N]),
  DnS = expression(Delta*bar(italic(n))[S]),
  V0 = expression(list(italic("V") * degree, "cm" ^ 3 ~ "mol" ^ -1)),
  DV0 = expression(list(Delta * italic("V") * degree, "cm" ^ 3 ~ "mol" ^ -1))
)
