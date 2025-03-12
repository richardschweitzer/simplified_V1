check_if_same_count_and_remove <- function(x, x_seed=NaN) {
  # get number of occurrences
  freq <- as.data.frame(table(x))
  colnames(freq) <- c("x", "Freq")
  # check whether the variance of different occurrences is not 0
  if (length(x)>1) {
    current_var = var(freq$Freq)
  } else {
    current_var = 0
  }
  # if the variance is not null, i.e., there are varying frequencies, continue
  remove_these <- NULL
  if (current_var!=0) {
    freq$Freq_diff <- freq$Freq - min(freq$Freq)
    # remove randomly chose values that have Freq_diff larger than 0
    for (r in 1:nrow(freq)) {
      freq_diff <- freq$Freq_diff[r]
      if (freq_diff > 0) {
        # add them to vector of items to remove
        if (!is.na(x_seed)) {
          set.seed(x_seed+r)
        }
        remove_these <- c(remove_these, sample(which(x==freq$x[r]))[1:freq_diff])
      }
    }
  }
  return(remove_these)
}