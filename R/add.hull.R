# canprot/add.hull.R
# Add convex hull around points 20200923
# Moved from chem16S::add_hull() 20240310

add.hull <- function(x, y = NULL, ...) {
  xy <- xy.coords(x, y, recycle = TRUE, setLab = FALSE)
  ihull <- chull(xy$x, xy$y)
  polygon(xy$x[ihull], xy$y[ihull], ...)
  invisible(ihull)
}
