# Text for figure labels
# Moved from internal.R and exported 20200204
# Moved from diffplot.R and changed expression() to quote() 20230617

cplab <- list(

  # Metrics in canprot
  Zc = quote(italic(Z)[C]),
  nH2O = quote(italic(n)[H[2]*O]),
  nO2 = quote(italic(n)[O[2]]),
  GRAVY = "GRAVY",
  pI = "pI",
  MW = "MW per residue",
  pMW = "MW per protein",
  V0 = quote("Volume per residue (cm" ^ 3 ~ "mol" ^ -1 * ")"),
  pV0 = quote("Volume per protein (cm" ^ 3 ~ "mol" ^ -1 * ")"),
  Density = "Density per residue",
  HC = "H/C",
  NC = "N/C",
  OC = "O/C",
  SC = "S/C",
  plength = "Protein length",
  Cost = "Metabolic cost",
  RespiratoryCost = "Respiratory cost",
  FermentativeCost = "Fermentative cost",
  B20Cost = "Biosynthetic cost in bacteria",
  Y20Cost = "Biosynthetic cost in yeast",
  H11Cost = "Biosynthetic cost in humans"

#  # Some others that aren't used (yet) in canprot
#  nAA = quote(italic(n)[AA]),
#  nC = quote(italic(n)[C] * "/AA"),
#  nH = quote(italic(n)[H] * "/AA"),
#  nN = quote(italic(n)[N] * "/AA"),
#  nO = quote(italic(n)[O] * "/AA"),
#  nS = quote(italic(n)[S] * "/AA"),
#  logfO2 = quote(log~italic("f")[O[2]]),
#  logaH2O = quote(log~italic("a")[H[2]*O])

)
