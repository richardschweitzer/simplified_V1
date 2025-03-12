compute_difference_3d_matrices <- function(mat_1, mat_2, 
                                           use_delayed_array="HDF5Array", 
                                           absolute_difference=TRUE) {
  # this function computes the absolute difference between two arrays [y,x,t]. 
  # Each may have different dimensions in [y,x], but they must be matched by dimnames.
  # In the end, we compute the mean for each time point. 
  require(assertthat)
  if (use_delayed_array == "RleArray") {
    require(DelayedArray)
  }
  if (use_delayed_array == "HDF5Array") {
    require(HDF5Array)
  }
  if (use_delayed_array == "SparseArray") {
    require(SparseArray)
  }
  
  # the easy case where both matrices have the same dimensions
  if (all(dim(mat_1)==dim(mat_2))) {
    
    # subtract them
    if (absolute_difference) {
      mat_diff <- abs(mat_2 - mat_1)
    } else {
      mat_diff <- mat_2 - mat_1
    }
    # perform mean-apply here
    res <- apply(X = mat_diff, MARGIN = 3, FUN = mean)
    rm(mat_diff)
    
  } else { # the more difficult case where the matrices have different dimensions
    
    # check whether matrices are the same in their 3rd (temporal dimension)
    mat_1_dimnames <- dimnames(mat_1)
    mat_2_dimnames <- dimnames(mat_2)
    are_equal(mat_1_dimnames[[3]], mat_2_dimnames[[3]])
    
    # check whether x-y dimnames are overlapping
    mat_1_dims <- dim(mat_1)
    mat_2_dims <- dim(mat_2)
    colnames_intersect <- intersect(mat_1_dimnames[[2]], mat_2_dimnames[[2]])
    rownames_intersect <- intersect(mat_1_dimnames[[1]], mat_2_dimnames[[1]])
    assert_that(length(colnames_intersect)>=10 & length(rownames_intersect)>=10, 
                msg = "There are less than 10 overlapping x-y dimnames in the two arrays!")
    # do they have the same resolution?
    are_equal(median(diff(as.numeric(mat_1_dimnames[[1]]))), median(diff(as.numeric(mat_2_dimnames[[1]]))))
    are_equal(median(diff(as.numeric(mat_1_dimnames[[2]]))), median(diff(as.numeric(mat_2_dimnames[[2]]))))
    spatial_resolution <- median(diff(as.numeric(mat_1_dimnames[[2]])))
    
    # create a large, padded array
    min_rowname <- min(c(as.numeric(mat_1_dimnames[[1]]), as.numeric(mat_2_dimnames[[1]])))
    max_rowname <- max(c(as.numeric(mat_1_dimnames[[1]]), as.numeric(mat_2_dimnames[[1]])))
    min_colname <- min(c(as.numeric(mat_1_dimnames[[2]]), as.numeric(mat_2_dimnames[[2]])))
    max_colname <- max(c(as.numeric(mat_1_dimnames[[2]]), as.numeric(mat_2_dimnames[[2]])))
    mat_size_cols <- seq(min_colname, max_colname, by = spatial_resolution)
    mat_size_rows <- seq(min_rowname, max_rowname, by = spatial_resolution)
    required_dims <- c(length(mat_size_rows), length(mat_size_cols), length(mat_1_dimnames[[3]]))
    required_dimnames <- list(mat_size_rows, mat_size_cols, mat_1_dimnames[[3]])
    if (use_delayed_array == "RleArray") {
      dat <- Rle(values = 0, lengths = prod(required_dims))
      padded_array <- RleArray(data = dat, 
                               dim = required_dims, dimnames = required_dimnames)
      rm(dat)
    } else if (use_delayed_array == "HDF5Array") { 
      # this is more memory-intense
      # padded_array <- DelayedArray(array(data = 0, 
      #                                    dim = required_dims, dimnames = required_dimnames))
      # this saves memory:
      padded_array <- as(array(data = 0,
                               dim = required_dims, dimnames = required_dimnames), "HDF5Array")
    } else if (use_delayed_array == "SparseArray") {
      padded_array <- SparseArray(x = array(data = 0, 
                                            dim = required_dims, dimnames = required_dimnames))
    } else {
      padded_array <- array(data = 0, 
                            dim = required_dims, dimnames = required_dimnames)
    }
    padded_dims <- dim(padded_array)
    padded_dimnames <- dimnames(padded_array)
    
    # enter mat_1
    mat_1_rowstart <- which.min(abs(as.numeric(padded_dimnames[[1]])-min(as.numeric(mat_1_dimnames[[1]]))))
    mat_1_colstart <- which.min(abs(as.numeric(padded_dimnames[[2]])-min(as.numeric(mat_1_dimnames[[2]]))))
    mat_1_padded <- padded_array
    if (use_delayed_array == "RleArray" | use_delayed_array == "HDF5Array" | use_delayed_array == "SparseArray") {
      mat_1_padded[mat_1_rowstart:(mat_1_rowstart+mat_1_dims[1]-1), 
                   mat_1_colstart:(mat_1_colstart+mat_1_dims[2]-1), ] <- mat_1
    } else {
      mat_1_padded[mat_1_rowstart:(mat_1_rowstart+mat_1_dims[1]-1), 
                   mat_1_colstart:(mat_1_colstart+mat_1_dims[2]-1), ] <- as.array(mat_1)
    }
    # enter mat_2
    mat_2_rowstart <- which.min(abs(as.numeric(padded_dimnames[[1]])-min(as.numeric(mat_2_dimnames[[1]]))))
    mat_2_colstart <- which.min(abs(as.numeric(padded_dimnames[[2]])-min(as.numeric(mat_2_dimnames[[2]]))))
    mat_2_padded <- padded_array
    if (use_delayed_array == "RleArray" | use_delayed_array == "HDF5Array" | use_delayed_array == "SparseArray") {
      mat_2_padded[mat_2_rowstart:(mat_2_rowstart+mat_2_dims[1]-1), 
                   mat_2_colstart:(mat_2_colstart+mat_2_dims[2]-1), ] <- mat_2
    } else {
      mat_2_padded[mat_2_rowstart:(mat_2_rowstart+mat_2_dims[1]-1), 
                   mat_2_colstart:(mat_2_colstart+mat_2_dims[2]-1), ] <- as.array(mat_2)
    }
    rm(padded_array)
    
    # checks
    are_equal(dim(mat_2_padded), dim(mat_1_padded))
    are_equal(dimnames(mat_2_padded), dimnames(mat_1_padded))
    # clean up
    rm(mat_1, mat_2)
    
    # compute difference for each time point
    if (absolute_difference) {
      mat_diff <- abs(mat_2_padded - mat_1_padded)
    } else {
      mat_diff <- mat_2_padded - mat_1_padded
    }
    # perform mean-apply here
    res <- apply(X = mat_diff, MARGIN = 3, FUN = mean)
    # clean up again
    rm(mat_1_padded, mat_2_padded, mat_diff)
  }
  
  return(res)
}

# # aux-function
# is.element.special <- function(x, y) {
#   if (length(y)==0) {
#     res <- rep(TRUE, times = length(x))
#   } else {
#     res <- is.element(x, y)
#   }
#   return(res)
# }



# plot_heatmap(as.matrix(mat_1[ , , 200]))
# plot_heatmap(as.matrix(mat_2[ , , 200]))
# plot_heatmap(as.matrix(mat_intersect_diff[ , , 200]))
# plot_heatmap(as.matrix(mat_1_only_diff[ , , 200]))
# plot_heatmap(as.matrix(mat_1_only_diff[ , , 400]))

# plot_heatmap(as.matrix(mat_1_padded[ , , 200]))
# plot_heatmap(as.matrix(mat_2_padded[ , , 200]))

# plot_heatmap(as.matrix(mat_1_padded[ , , 500]))
# plot_heatmap(as.matrix(mat_2_padded[ , , 500]))



