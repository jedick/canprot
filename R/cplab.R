# Text for figure labels
# Moved from internal.R and exported 20200204
# Moved from diffplot.R and changed expression() to quote() 20230617

cplab <- list(
  nH2O = quote(italic(n)[H[2]*O]),
  DnH2O = quote(Delta*italic(n)[H[2]*O]),
  nO2 = quote(italic(n)[O[2]]),
  DnO2 = quote(Delta*italic(n)[O[2]]),
  Zc = quote(italic(Z)[C]),
  DZc = quote(Delta*italic(Z)[C]),
  nC = quote(italic(n)[C] * "/AA"),
  nN = quote(italic(n)[N] * "/AA"),
  nS = quote(italic(n)[S] * "/AA"),
  DnC = quote(Delta*italic(n)[C] * "/AA"),
  DnN = quote(Delta*italic(n)[N] * "/AA"),
  DnS = quote(Delta*italic(n)[S] * "/AA"),
  V0 = quote("Volume per residue (cm" ^ 3 ~ "mol" ^ -1 * ")"),
  pV0 = quote("Volume per protein (cm" ^ 3 ~ "mol" ^ -1 * ")"),
  DV0 = quote(list(Delta * italic("V") * degree, "cm" ^ 3 ~ "mol" ^ -1)),
  nAA = quote(italic(n)[AA]),
  DnAA = quote(Delta*italic(n)[AA]),
  GRAVY = "GRAVY",
  DGRAVY = quote(Delta*"GRAVY"),
  pI = "pI",
  DpI = quote(Delta*"pI"),
  MW = quote("MW per residue"),
  DMW = quote(Delta * "MW"),
  plength = "Protein length"
)
