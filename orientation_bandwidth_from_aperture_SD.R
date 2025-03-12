orientation_bandwidth_from_aperture_SD <- function(SF, Gaussian_SD) {
  # determine the orientation bandwidth of a Gabor filter based on its SF and the aperture SD.
  # this is all according to Movellan "Tutorial on Gabor Filters"
  # https://inc.ucsd.edu/mplab/75/media//gabor.pdf
  C <- sqrt(log(2)/pi)
  b <- sqrt(1/(2*pi*Gaussian_SD^2))
  ori_band <- 2 * atan((b*C)/SF) # this is the full bandwidth, not half!
  return(ori_band)
}