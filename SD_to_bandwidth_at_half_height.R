SD_to_bandwidth_at_half_height <- function(SD, run_inverse = FALSE) { 
  # see: https://en.wikipedia.org/wiki/Full_width_at_half_maximum
  if (!run_inverse) {
    FWHM = 2 * sqrt(2 * log(2)) * SD
  } else {
    FWHM = SD / (2 * sqrt(2 * log(2))) # here SD is FWHM and vice verse
    # so, actually: SD = FWHM / (2 * sqrt(2 * log(2)))
  }
  return(FWHM)
}