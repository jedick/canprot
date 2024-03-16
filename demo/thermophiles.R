# canprot/demo/thermophiles.R
# Distinguishing thermophile genomes using entropy and oxidation state
# 20240310

library(canprot)

## @knitr thermophiles_demo_body

# Change this to TRUE to make the image for the README
# (with only the Nitrososphaeria plot)
do.pdf <- FALSE
if(do.pdf) pdf("thermophiles.pdf", width = 4.2, height = 4.2)
if(!do.pdf) par(mfrow = c(1, 2))

par(mar = c(4, 4.5, 2, 1))

# Define metrics to calculate
metrics <- c("Zc", "S0g")
xlab <- cplab[[metrics[1]]]
ylab <- cplab[[metrics[2]]]

# Define common y-axis limit
ylim <- c(0.3000, 0.3045)

# Define colors
col2.8 <- adjustcolor(2, alpha.f = 0.8)
col4.8 <- adjustcolor(4, alpha.f = 0.8)
col2.5 <- adjustcolor(2, alpha.f = 0.5)

if(!do.pdf) {

## Plot 1: methanogen genomes

# Read amino acid composition
aafile <- system.file("extdata/aa/methanogen_aa.csv", package = "canprot")
aa <- read.csv(aafile)
# Set fill color for thermophiles
ithermo <- aa$ref > 50
bg <- ifelse(ithermo, col2.8, col4.8)
# Set point symbol for Class I and II
pch <- ifelse(aa$abbrv == "Class I", 21, 22)
# Calculate and plot metrics
metvals <- calc_metrics(aa, metrics)
plot(metvals, pch = pch, bg = bg, xlab = xlab, ylab = ylab, ylim = ylim)
# Add convex hull around thermophiles
add_hull(metvals[ithermo, ], col = col2.5, border = NA)
# Add legend and title
text(-0.225, 0.303, "Thermophiles", col = 2)
text(-0.18, 0.3003, "Mesophiles", col = 4)
legend("bottomleft", c("Class I", "Class II"), pch = c(21, 22), cex = 0.9)
title("Methanogen genomes", font.main = 1)

}

## Plot 2: Nitrososphaeria MAGs

# Read genome data
datfile <- system.file("extdata/aa/nitrososphaeria_MAGs.csv", package = "canprot")
dat <- read.csv(datfile)
# Read amino acid compositions
aafile <- system.file("extdata/aa/nitrososphaeria_aa.csv", package = "canprot")
aa <- read.csv(aafile)
# Match accessions
idat <- match(aa$organism, dat$Accession)
dat <- dat[idat, ]
# Filter by family
ifam <- !dat$Family %in% c("Nitrosopumilaceae", "Nitrososphaeraceae")
dat <- dat[ifam, ]
aa <- aa[ifam, ]
# Assign point symbol and color
pch <- sapply(dat$Respiration.type, switch, Anaerobic = 21, Aerobic = 22, Microaerobic = 23)
bg <- sapply(dat$Habitat.type, switch, Thermal = col2.8, Nonthermal = col4.8)
# Calculate and plot metrics
metvals <- calc_metrics(aa, metrics)
# Make plot
plot(metvals, pch = pch, bg = bg, xlab = xlab, ylab = ylab, ylim = ylim)
# Add convex hull around MAGs from thermal habitats
ithermal <- dat$Habitat.type == "Thermal"
add_hull(metvals[ithermal, ], col = col2.5, border = NA)
# Add legend and title
text(-0.215, 0.303, "Thermal habitats", col = 2)
text(-0.175, 0.3013, "Nonthermal habitats", col = 4)
legend("bottomright", c("Anaerobic", "Aerobic", "Microaerobic"), pch = c(21, 22, 23), cex = 0.9)
title("Nitrososphaeria MAGs", font.main = 1)

if(do.pdf) dev.off()
