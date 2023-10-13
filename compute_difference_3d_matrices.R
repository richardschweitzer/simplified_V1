compute_difference_3d_matrices <- function(mat_1, mat_2, 
                                           use_RLE=TRUE, 
                                           absolute_difference=TRUE) {
  # this function computes the absolute difference between two arrays [y,x,t]. 
  # Each may have different dimensions in [y,x], but they must be matched by dimnames.
  # In the end, we compute the mean for each time point. 
  require(assertthat)
  if (use_RLE) {
    require(DelayedArray)
  }
  
  # check whether matrices are the same in their 3rd (temporal dimension)
  mat_1_dimnames <- dimnames(mat_1)
  mat_2_dimnames <- dimnames(mat_2)
  are_equal(mat_1_dimnames[[3]], mat_2_dimnames[[3]])
  
  # check whether x-y dimnames are overlapping
  mat_1_dims <- dim(mat_1)
  mat_2_dims <- dim(mat_2)
  colnames_intersect <- intersect(mat_1_dimnames[[2]], mat_2_dimnames[[2]])
  rownames_intersect <- intersect(mat_1_dimnames[[1]], mat_2_dimnames[[1]])
  assert_that(length(colnames_intersect)>10 & length(rownames_intersect)>10, 
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
  if (use_RLE) {
    dat <- Rle(values = 0, lengths = prod(required_dims))
    padded_array <- RleArray(data = dat, 
                             dim = required_dims, dimnames = required_dimnames)
    rm(dat)
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
  if (use_RLE) {
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
  if (use_RLE) {
    mat_2_padded[mat_2_rowstart:(mat_2_rowstart+mat_2_dims[1]-1), 
                 mat_2_colstart:(mat_2_colstart+mat_2_dims[2]-1), ] <- mat_2
  } else {
    mat_2_padded[mat_2_rowstart:(mat_2_rowstart+mat_2_dims[1]-1), 
                 mat_2_colstart:(mat_2_colstart+mat_2_dims[2]-1), ] <- as.array(mat_2)
  }
  
  # checks
  are_equal(dim(mat_2_padded), dim(mat_1_padded))
  are_equal(dimnames(mat_2_padded), dimnames(mat_1_padded))
  
  # compute difference for each time point
  if (absolute_difference) {
    res <- apply(X = abs(mat_2_padded - mat_1_padded), 
                 MARGIN = 3, FUN = mean)
  } else {
    res <- apply(X = mat_2_padded - mat_1_padded, 
                 MARGIN = 3, FUN = mean)
  }
  
  # else { # in this version, we split the arrays up, which should need less memory
  #   
  #   mat_1_rownames_only <- setdiff(mat_1_dimnames[[1]], mat_2_dimnames[[1]]) # only mat_1 has those dimnames
  #   mat_1_colnames_only <- setdiff(mat_1_dimnames[[2]], mat_2_dimnames[[2]])
  #   mat_2_rownames_only <- setdiff(mat_2_dimnames[[1]], mat_1_dimnames[[1]]) # only mat_2 has those dimnames
  #   mat_2_colnames_only <- setdiff(mat_2_dimnames[[2]], mat_1_dimnames[[2]])
  #   
  #   # first, compute difference for matching dimnames
  #   mat_intersect_diff <- mat_1[is.element(mat_1_dimnames[[1]], rownames_intersect), 
  #                               is.element(mat_1_dimnames[[2]], colnames_intersect), , drop = FALSE] - 
  #     mat_2[is.element(mat_2_dimnames[[1]], rownames_intersect), 
  #           is.element(mat_2_dimnames[[2]], colnames_intersect), , drop = FALSE]
  #   if (absolute_difference) {
  #     mat_intersect_diff <- abs(mat_intersect_diff)
  #   }
  #   
  #   ## second, compute difference for non-matching dimnames
  #   mat_1_only_diff <- mat_1[is.element.special(mat_1_dimnames[[1]], mat_1_rownames_only), 
  #                            is.element.special(mat_1_dimnames[[2]], mat_1_colnames_only), , drop = FALSE] - 0
  #   if (absolute_difference) {
  #     mat_1_only_diff <- abs(mat_1_only_diff)
  #   }
  #   mat_2_only_diff <- mat_2[is.element.special(mat_2_dimnames[[1]], mat_2_rownames_only), 
  #                            is.element.special(mat_2_dimnames[[2]], mat_2_colnames_only), , drop = FALSE] - 0
  #   if (absolute_difference) {
  #     mat_2_only_diff <- abs(mat_2_only_diff)
  #   }
  #   
  #   # compute difference for each time point
  #   if (absolute_difference) {
  #     res <- apply(X = mat_1_only_diff + mat_1_only_diff + mat_2_only_diff, 
  #                  MARGIN = 3, FUN = mean)
  #   } else {
  #     res <- apply(X = mat_2_padded - mat_1_padded, 
  #                  MARGIN = 3, FUN = mean)
  #   }
  #   
  # }
  # 
  # 
  
  return(res)
}

# aux-function
is.element.special <- function(x, y) {
  if (length(y)==0) {
    res <- rep(TRUE, times = length(x))
  } else {
    res <- is.element(x, y)
  }
  return(res)
}



# plot_heatmap(as.matrix(mat_1[ , , 200]))
# plot_heatmap(as.matrix(mat_2[ , , 200]))
# plot_heatmap(as.matrix(mat_intersect_diff[ , , 200]))
# plot_heatmap(as.matrix(mat_1_only_diff[ , , 200]))
# plot_heatmap(as.matrix(mat_1_only_diff[ , , 400]))

# plot_heatmap(as.matrix(mat_1_padded[ , , 200]))
# plot_heatmap(as.matrix(mat_2_padded[ , , 200]))

# plot_heatmap(as.matrix(mat_1_padded[ , , 500]))
# plot_heatmap(as.matrix(mat_2_padded[ , , 500]))



