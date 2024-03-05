# canprot/demo/redoxin.R
# Zc of ferredoxin, thioredoxin, and glutaredoxin vs midpoint reduction potential
# Based on Fig. 5 in this preprint: http://dx.doi.org/10.1101/004127 dated 20140414
# Moved from JMDplots::aoscp99() to canprot on 20240304
#   and modified to remove PDI (protein disulfide isomerase),
#   which was incorrectly grouped with proteins from E. coli

library(canprot)

## @knitr redoxin_demo_body

# Data file with protein IDs, sequence start/stop positions, and midpoint potentials
data_file <- system.file("extdata/fasta/redoxin.csv", package = "canprot")
dat <- read.csv(data_file)
# Drop PDI (a human protein) 20240304
dat <- dat[dat$protein != "PDI", ]

# Read header lines
fasta_file <- system.file("extdata/fasta/redoxin.fasta", package = "canprot")
# TODO: change canprot::read.fasta() to read.fasta() after CHNOSZ 2.1.0 is superseded on CRAN 20240304
headers <- canprot::read.fasta(fasta_file, type = "header")
# Locate the sequences in the FASTA file
iseqs <- sapply(dat$ID, grep, x = headers)
# Loop over proteins
aalist <- lapply(1:nrow(dat), function(i) {
  # Read the amino acid composition of this protein
  canprot::read.fasta(fasta_file, iseq = iseqs[i], start = dat$start[i], stop = dat$stop[i])
})
aa <- do.call(rbind, aalist)

# Make ferredoxin-thioredoxin reductase dimer (variable chain/catalytic chain)
iFTR <- grep("FTR", dat$protein)
aa[iFTR[1], 6:24] <- colSums(aa[iFTR, 6:24])
aa <- aa[-iFTR[2], ]
dat$protein[iFTR[1]] <- paste(dat$protein[iFTR], collapse = ":")
dat <- dat[-iFTR[2], ]

# Calculate Zc
Zc_values <- Zc(aa)

# Point symbols for E. coli and spinach
pch <- rep(19, length(Zc_values))
pch[dat$organism=="spinach"] <- 0

# Start plot
par(las = 1)
plot(dat$E0, Zc_values, pch = pch, xlim = c(-450, -100), ylim = c(-0.28, -0.04),
  xlab = expression(list(italic(E)*degree*"'", mV)),
  ylab = expression(italic(Z)[C]))
# Add dashed lines
lines(dat$E0[dat$organism == "ecoli"], Zc_values[dat$organism == "ecoli"], lty = 2)
lines(dat$E0[dat$organism == "spinach"], Zc_values[dat$organism == "spinach"], lty = 2)

# Add labels
pos <- rep(1, length(Zc_values))
pos[dat$organism == "ecoli"] <- 4
pos[dat$organism == "spinach"] <- 2
dx <- dy <- numeric(length(Zc_values))
dx[dat$protein == "DsbA"] <- -10
dy[dat$protein == "DsbA"] <- -0.012
dx[dat$protein == "DsbC"] <- -20
dy[dat$protein == "DsbC"] <- -0.012
text(dat$E0 + dx, Zc_values + dy, dat$protein, pos = pos)

# Add legend
legend("bottomright", pch = c(19, 0), legend = c(expression(italic("E. coli")), "spinach"))
