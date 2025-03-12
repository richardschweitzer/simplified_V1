make_eeg_like_downsample <- function(mat_2d,    # 2D matrix
                                     integrate_sd = 75, # in pixels, SD of the integration window
                                     n_points_xy = c(6, 10), # number of sample points in vertical and horizontal dimensions
                                     d_points_xy = c(40, 40),  # distance [pixels] between sample points in vertical and horizontal dimensions 
                                     do_plot = FALSE
) {
  require(assertthat)
  require(pracma)
  require(fields) # for interpolation
  # compute the grids
  x_points = c(rev(-1 * seq(0, d_points_xy[2]*n_points_xy[2]/2, length.out = n_points_xy[2]/2)[-1]), 
                seq(0, d_points_xy[2]*n_points_xy[2]/2, length.out = n_points_xy[2]/2))
  y_points = c(rev(-1 * seq(0, d_points_xy[1]*n_points_xy[1]/2, length.out = n_points_xy[1]/2)[-1]), 
               seq(0, d_points_xy[1]*n_points_xy[1]/2, length.out = n_points_xy[1]/2))
  xy_grid <- expand.grid(list(x = x_points, y = y_points))
  xy_grid$sensor_name <- paste(round(xy_grid$x, 1), round(xy_grid$y, 1), sep = "_")
  # make the integration window
  # create the meshgrid and the Gaussian integration window
  integrate_seq <- c(rev(-1*seq(0, 3*integrate_sd, 1)[-1]), seq(0, 3*integrate_sd, 1))
  meshlist <- meshgrid(x = integrate_seq, y = integrate_seq) 
  gaussian_field <- exp(-(meshlist$X^2+meshlist$Y^2)/(2*integrate_sd^2))
  # make sure mat_2d has proper dimnames
  assert_that(!is.null(dimnames(mat_2d)), msg = "mat_2d must have proper dimnames specifying the pixel dimensions!")
  x_res = as.numeric(dimnames(mat_2d)[[2]])
  y_res = as.numeric(dimnames(mat_2d)[[1]])
  assert_that(all(!is.na(c(x_res, y_res))), msg = "dimnames of mat_2d contain values not convertable to numeric!")
  # then make the interpolation object
  image_obj <- list(x = y_res, y = x_res, z = mat_2d) 
  ## now go through points
  xy_grid$output <- NaN
  for (xy_i in 1:nrow(xy_grid)) {
    x_now = xy_grid$x[xy_i]
    y_now = xy_grid$y[xy_i]
    # which points to select?
    select_col <- x_now + integrate_seq
    select_row <- y_now + integrate_seq
    query_points <- expand.grid(rows = select_row, cols = select_col)
    # select current subimage via linear interpolation (requires package fields)
    target_image <- interp.surface( obj = image_obj, loc = query_points )
    target_image <- matrix(data = target_image, nrow = length(select_row), ncol = length(select_col),
                           byrow = FALSE, dimnames = list(select_row, select_col))
    target_image[is.na(target_image)] <- 0
    # to check: plot_heatmap(gaussian_field)
    # third, compute weighted average
    wa <- weighted.mean(x = target_image, w = gaussian_field) * 1000
    # finally, save output
    xy_grid$output[xy_i] <- wa
  }
  # plot?
  if (do_plot) {
    require(ggplot2)
    require(viridis)
    mat_long <- reshape2::melt(mat_2d)
    colnames(mat_long) <- c("y", "x", "z")
    p <- ggplot(mat_long, aes(x = x, y = y)) + 
      geom_tile(aes(fill = z)) + 
      scale_y_reverse() + theme_minimal() + 
      coord_cartesian(expand = FALSE) + 
      geom_point(data = xy_grid, aes(x = x, y = y, color = output)) + 
      scale_fill_viridis_c(option = "mako") + 
      scale_color_viridis_c(option = "magma", begin = 0.2)
    print(p)
  }
  # return output
  return(xy_grid)
}



