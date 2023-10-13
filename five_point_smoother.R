five_point_smoother <- function(x) { # 5-point moving window
  N <- length(x)
  v <- vector(mode = "numeric", length = N)
  for (r in 1:N) {
    if (r==N || r == 1) { 
      v[r] = x[r]
    } else if (r==N-1 || r == 2) {
      v[r] = (x[r-1] + x[r] + x[r+1]) / 3
    } else {
      v[r] = (x[r-2] + x[r-1] + x[r] + x[r+1] + x[r+2]) / 5
    }
  }
  return(v)
}