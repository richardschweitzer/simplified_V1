pad_and_extrapolate <- function(x, smooth_iterations=3, average_filter_size=5, do_plot=FALSE) {
  require(smoothie)
  # standard params:
  # smooth_iterations <- 3 # perform smooth and reinsert how many times?
  # average_filter_size <- 5 # kernel dimensions of average filter
  
  # get noise patch
  noise_temp <- x
  current_size <- dim(noise_temp)
  # now compute pad factor based on the parameters above
  pad_factor <- (current_size + smooth_iterations * average_filter_size) / current_size
  future_size <- ceiling(current_size * pad_factor)
  # check whether x has dimnames
  if (!is.null(dimnames(noise_temp))) {
    require(Hmisc)
    require(assertthat)
    has_dimnames <- TRUE
    # extrapolate dimnames
    future_dimnames <- vector(mode = "list", length = 2)
    for (dim_i in (1:2)) {
      current_dimnames_1 <- as.numeric(dimnames(noise_temp)[[dim_i]])
      current_indeces_1 <- ceiling(future_size[dim_i]/2-current_size[dim_i]/2):
        (ceiling(future_size[dim_i]/2-current_size[dim_i]/2)+current_size[dim_i]-1)
      assert_that(length(current_dimnames_1)==length(current_indeces_1))
      future_indeces_1 <- 1:future_size[dim_i]
      future_dimnames_1 <- approxExtrap(current_indeces_1, current_dimnames_1, future_indeces_1)
      assert_that(length(future_dimnames_1$y)==length(future_indeces_1))
      future_dimnames[[dim_i]] <- future_dimnames_1$y
    }
  } else {
    has_dimnames <- FALSE
  }
  # pad the noise
  future_noise <- matrix(0, ncol = future_size[2], nrow = future_size[1])
  future_noise[ceiling(future_size[1]/2-current_size[1]/2):(ceiling(future_size[1]/2-current_size[1]/2)+current_size[1]-1), 
               ceiling(future_size[2]/2-current_size[2]/2):(ceiling(future_size[2]/2-current_size[2]/2)+current_size[2]-1)] <-
    noise_temp
  if (do_plot) {
    par(mfrow=c(1,smooth_iterations+1))
    plot(diff(future_noise[future_size[1]/2, ]), type = "l", main = "original edges") # here you nicely see the edges
  }
  for (smooth_iteration in 1:smooth_iterations) {
    # smooth the noise
    future_noise <- kernel2dsmooth(x = future_noise, kernel.type = "average", 
                                   nx=average_filter_size, ny=average_filter_size)
    # reinsert original noise
    future_noise[ceiling(future_size[1]/2-current_size[1]/2):(ceiling(future_size[1]/2-current_size[1]/2)+current_size[1]-1), 
                 ceiling(future_size[2]/2-current_size[2]/2):(ceiling(future_size[2]/2-current_size[2]/2)+current_size[2]-1)] <-
      noise_temp
    #plot_heatmap(future_noise)
    if (do_plot) {
      plot(diff(future_noise[future_size[1]/2, ]), type = "l", main = paste("edge adjustment", smooth_iteration))
    }
  }
  # apply dimnames?
  if (has_dimnames) {
    dimnames(future_noise) <- future_dimnames
  }
  return(future_noise)
}