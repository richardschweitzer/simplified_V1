bootstrap_sum_se <- function(x, n_reps = 1000, n_cores = 1, print_output = FALSE) {
  if (n_cores > 1) { 
    parallelize <- "multicore" 
  } else {
    parallelize <- "no"    
  }
  # compute summary
  med = sum(x, na.rm = TRUE) # set aggregation function: median, mean, sum, whatever
  if (length(x)>0 & !is.na(med)) {
    # for bootstrapping we need a sampling function
    samplemean <- function(x, d) {
      return(sum(x[d], na.rm = TRUE)) # set aggregation function: median, mean, sum, whatever
    }
    # now we bootstrap to estimate the median
    require(boot)
    booty = boot(x, samplemean, 
                 R = n_reps, parallel = parallelize, ncpus = n_cores)  
    # compute 
    se <- sd(booty$t)
    # compute confidence interval based on Efron's confident limit
    upper <- quantile(booty$t, 0.975, names = FALSE)
    lower <- quantile(booty$t, 0.025, names = FALSE)
  } else {
    se <- NA
    upper <- NA
    lower <- NA
  }
  # return values
  res <- list(mean = med, se = se, lower = lower, upper = upper)
  if (print_output) {
    print(res)
  }
  return(res)
}
