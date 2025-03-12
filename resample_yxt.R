resample_yxt <- function(a_3d, target_size_yxt, 
                         interpolate_method = "trilinear", # available are: nearest, trilinear, area
                         preserve_corners = TRUE) {
  # input should be a [y,x,t] torch float tensor that should be resampled to target size [y2,x2,t2]
  require(torch)
  require(assertthat)
  # make sure constraints are met
  assert_that(length(dim(a_3d))==3 & length(target_size_yxt)==3)
  # unsqueeze to make it 5D
  a_3d$unsqueeze_(1)
  a_3d$unsqueeze_(1)
  # perform interpolation here
  a_3d <- nnf_interpolate(input = a_3d, 
                          size = target_size_yxt, 
                          mode = interpolate_method, 
                          align_corners = preserve_corners)
  # squeeze again
  a_3d$squeeze_(1)
  a_3d$squeeze_(1)
  # done
  return(a_3d)
}