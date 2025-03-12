get_gabor_field_heiko <- function(rf_freq_dva, rf_ori, 
                                  bw = c(0.5945, 0.2965), # according to Heiko's paper
                                  im_size_dva, 
                                  scr.ppd, ppd_scaler = 1 # the ppd scaler can be applied to produce more highres gabor filters
                                  ) {
  # uses Heiko Schuett's function to create zero-phase and half-phase Gabor filters
  im_size_dva <- c(im_size_dva, im_size_dva)
  im_size_pix <- round(im_size_dva*scr.ppd)
  # check for odd dimensions
  if (!mod(im_size_pix[1], 2) == 1) {
    im_size_pix[1] <- im_size_pix[1] + 1
    im_size_dva[1] <- im_size_pix[1] / scr.ppd
  }
  if (!mod(im_size_pix[2], 2) == 1) {
    im_size_pix[2] <- im_size_pix[2] + 1
    im_size_dva[2] <- im_size_pix[2] / scr.ppd
  }
  # get those
  rf_field_zero <- get_log_gabor_heiko(imSize = im_size_pix, 
                                       degSize = im_size_dva, 
                                       bw = bw, # tuple for SF and Ori
                                       freq = rf_freq_dva, orientation = rf_ori, phase = 0, 
                                       scaler = ppd_scaler)
  rf_field_half <- get_log_gabor_heiko(imSize = im_size_pix, 
                                       degSize = im_size_dva, 
                                       bw = bw,
                                       freq = rf_freq_dva, orientation = rf_ori, phase = pi*1/2, 
                                       scaler = ppd_scaler)
  # return data
  assert_that(all(dim(rf_field_half)==dim(rf_field_zero)) & dim(rf_field_zero)[1]==dim(rf_field_zero)[2] )
  rf_info <- c(rf_freq_dva = rf_freq_dva, 
               rf_mesh_width = round(im_size_dva*scr.ppd), 
               rf_ori = rf_ori, 
               scr.ppd = scr.ppd, ppd_scaler = ppd_scaler)
  return(list(rf_field_zero, rf_field_half, rf_info))
}