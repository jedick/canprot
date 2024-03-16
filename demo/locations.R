# canprot/demo/locations.R
# Plot chemical metrics of human proteins in different subcellular locations
# 20240302

library(canprot)

## @knitr locations_demo_setup

# Change to TRUE to make a PDF image
do.pdf <- FALSE
if(do.pdf) pdf("locations.pdf", width = 8, height = 4.5)

## @knitr locations_demo_body

# Read SI table
file <- system.file("extdata/protein/TAW+17_Table_S6_Validated.csv", package = "canprot")
dat <- read.csv(file)
# Keep only proteins with validated location
dat <- dat[dat$Reliability == "Validated", ]
# Keep only proteins with one annotated location
dat <- dat[rowSums(dat[, 4:32]) == 1, ]

# Get the amino acid compositions
aa <- human_aa(dat$Uniprot)
# Put the location into the amino acid data frame
aa$location <- dat$IF.main.protein.location

# Use top locations (and their colors) from Fig. 2B of Thul et al., 2017
locations <- c("Cytosol","Mitochondria","Nucleoplasm","Nucleus","Vesicles","Plasma membrane")
col <- c("#194964", "#2e6786", "#8a2729", "#b2333d", "#e0ce1d", "#e4d71c")
# Keep the proteins in these locations
aa <- aa[aa$location %in% locations, ]
## Keep only proteins with length between 100 and 2000
#aa <- aa[plength(aa) >= 100 & plength(aa) <= 2000, ]

# Get amino acid composition for proteins in each location
# (Loop over groups by piping location names into lapply)
aalist <- lapply(locations, function(location) aa[aa$location == location, ] )

# Setup plot
par(mfrow = c(1, 2))
titles <- c(Zc = "Carbon oxidation state", pI = "Isoelectric point")
# Calculate Zc and pI
for(metric in c("Zc", "pI")) {
  datlist <- lapply(aalist, metric)
  bp <- boxplot(datlist, ylab = cplab[[metric]], col = col, show.names = FALSE)
  add_cld(datlist, bp)
  # Make rotated labels
  x <- (1:6) + 0.1
  y <- par()$usr[3] - 1.5 * strheight("A")
  text(x, y, locations, srt = 25, adj = 1, xpd = NA)
  axis(1, labels = FALSE)
  title(titles[metric], font.main = 1)
}

## @knitr locations_demo_wrapup

if(do.pdf) dev.off()

# The PDF file in inst/images was converted to PNG with GIMP and compressed with pngquant
