# canprot/R/internal.R
# internal (non-exported) objects
# 20160706 jmd

# setup basis species
setbasis <- function(basis="QEC") {
  if(basis=="QEC+") {
    # QEC with H+
    basis(c("glutamine", "glutamic acid", "cysteine", "H2O", "oxygen", "H+"))
    basis("C5H10N2O3", -3.2)
    basis("C5H9NO4", -4.5)
    basis("C3H7NO2S", -3.6)
    basis("H+", -7)
  } else if(basis=="QEC") {
    # amino acids (approximate concentrations in human plasma)
    basis(c("glutamine", "glutamic acid", "cysteine", "H2O", "oxygen"))
    basis("C5H10N2O3", -3.2)
    basis("C5H9NO4", -4.5)
    basis("C3H7NO2S", -3.6)
  } else if(basis=="QEC4") {
    # amino acids (basis set II from Dick, 2016)
    basis(c("glutamine", "glutamic acid", "cysteine", "H2O", "oxygen"))
    basis("C5H10N2O3", -4)
    basis("C5H9NO4", -4)
    basis("C3H7NO2S", -4)
  } else if(basis=="inorganic") {
    # inorganic species (basis set I from Dick, 2016)
    basis(c("CO2", "H2O", "NH3", "H2S", "oxygen"))
    basis(c("CO2", "NH3", "H2S"), c(-3, -4, -7))
  } else stop(paste("undefined basis setting:", basis))
}

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
  DnN = expression(Delta*bar(italic(n))[N]),
  DnS = expression(Delta*bar(italic(n))[S])
)
