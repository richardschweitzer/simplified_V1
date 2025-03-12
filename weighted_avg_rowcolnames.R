weighted_avg_rowcolnames <- function(x) {
  # x must be a matrix with nonnegative weights
  require(torch)
  a_rownames <- torch_tensor(matrix(data = as.numeric(rownames(x)), 
                                    nrow = length(rownames(x))))
  a_colnames <- torch_tensor(matrix(data = as.numeric(colnames(x)), 
                                    ncol = length(colnames(x))) )
  a <- torch_tensor(matrix(data = x, nrow = nrow(x), ncol = ncol(x), 
                           dimnames = list(as.numeric(a_rownames), 
                                           as.numeric(a_colnames)) ) )
  # use multiply that uses broadcasting, then sum over matrix
  weighted_avg_row <- as.numeric(
    torch_sum(torch_multiply(a, a_rownames), dim = c(1,2)) / 
      torch_sum(a, dim = c(1,2)) )
  weighted_avg_col <- as.numeric(
    torch_sum(torch_multiply(a, a_colnames), dim = c(1,2)) / 
      torch_sum(a, dim = c(1,2)) )
  rm(a, a_colnames, a_rownames)
  return(list(weighted_avg_row, weighted_avg_col))
}