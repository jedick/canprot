# canprot/R/get_colors.R
# get colors for rank-difference (potential) diagrams
# 20160710 jmd

get_colors <- function(x, max50=FALSE) {
  # diverging (blue - light grey - red) palette
  # max50: values over 50% are all deepest color (red or blue)
  if(max50) dcol <- colorspace::diverge_hcl(1000, c = 100, l = c(50, 90), power = 1)
  else dcol <- colorspace::diverge_hcl(2000, c = 100, l = c(50, 90), power = 1)
  # the range of values
  xrange <- range(x)
  # select range of colors corresponding to values
  if(any(xrange < 0)) {
    # white to deep blue
    if(max50) blues <- rev(c(rep(dcol[1], 500), dcol[1:500]))
    else blues <- rev(dcol[1:1000])
    blues <- blues[1:-round(10*xrange[1])]
  } else blues <- character()
  if(any(xrange > 0)) {
    # white to deep red
    if(max50) reds <- c(dcol[501:1000], rep(dcol[1000], 500))
    else reds <- dcol[1001:2000]
    reds <- reds[1:round(10*xrange[2])]
  } else reds <- character()
  col <- c(rev(blues), reds)
  return(col)
}
