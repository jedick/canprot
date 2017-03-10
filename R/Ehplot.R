# canprot/R/Ehplot.R
# show values of Eh as a function of logaH2O/logfO2
# 20160710 jmd

Ehplot <- function(basis="QEC", T=37, pH=7.4) {
  logfO2 <- c(-75, -55)
  H2O <- c(-10, 10)
  # logK for the reaction H2O(liq) = 2H+ + 2e- + 0.5O2(g)
  logK <- subcrt(c("H2O", "H+", "e-", "oxygen"), c(-1, 2, 2, 0.5), T=T)$out$logK
  # to calculate logaH2O at a given logfO2 and Eh
  logaH2O <- function(logfO2, Eh) {
    pe <- convert(Eh, "pe", T=convert(T, "K"))
    return(0.5*logfO2 - 2*pH - 2*pe - logK)
  }
  plot(0, 0, xlim=logfO2, ylim=H2O, xlab=cplab$logfO2, ylab=cplab$logaH2O, type="n", xaxs="i", yaxs="i")
  for(Eh in seq(-0.8, 0.2, by=0.2)) {
    lines(logfO2, logaH2O(logfO2, Eh)) 
    if(basis=="QEC") text(-61+19*Eh, logaH2O(-61+19*Eh, Eh) + 1, Eh)
    if(basis=="inorganic") text(-71+19*Eh, logaH2O(-71+19*Eh, Eh) + 1, Eh)
  }
  title(main="Eh (volt)", cex.main=1.1)
}
