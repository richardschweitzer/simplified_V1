anderson_burr_RF_width <- function(SF) {
  widths <- rep(NaN, length(SF))
  widths[SF<1] <- exp(-0.5 * log(SF[SF<1]) + log(1))
  widths[SF>=1] <- exp(-1 * log(SF[SF>=1]) + log(1)) # or simply: 1/SF[SF>=1]
  return(widths)
}