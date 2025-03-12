baseline_correction <- function(EEG_subset, columns_used, id_column, time_column, baseline_interval = c(-100, 0), 
                                    use_data_table = FALSE) {
  # baseline correction for SDEC
  # by richard, 09/2022
  
  # to debug:
  # id_column = "ID"
  # columns_used = signal_cols
  # time_column = "time_stim_on_bins"
  # baseline_interval = c(-100, 0)
  
  require(data.table)
  #require(parallel) # for mcl
  
  # should we use the data table or data frame version
  if (use_data_table) {
    warning('using baseline_correction_par with use_data_table=TRUE! This mode does not know how to deal with NaN in baseline period!')
    setDT(EEG_subset)
  } else {
    setDF(EEG_subset)  # convert to data.frame
  }
  
  # first, split data table into list
  if (use_data_table) {
    EEG_subset_list <- split(x = EEG_subset, by = eval(id_column))
  } else {
    EEG_subset_list <- split.data.frame(x = EEG_subset, f = EEG_subset[ , eval(id_column)])
  }
  rm(EEG_subset)
  
  # now, perform a mclapply (not strictly needed) on that list
  EEG_subset_list <- lapply(X = EEG_subset_list, 
                            FUN = function(x, columns_used, time_column, baseline_interval, use_data_table) {
                              # find out where the baseline is
                              if (use_data_table) {
                                require(data.table)
                                # get all the baseline means
                                baseline_time_here <- as.vector(t(x[ , ..time_column])) >= baseline_interval[1] & 
                                  as.vector(t(x[ , ..time_column])) <= baseline_interval[2]
                                baseline_means <- x[baseline_time_here, lapply(.SD, mean), .SDcols = columns_used]
                              } else {
                                baseline_time_here <- as.vector(t(x[ , eval(time_column)])) >= baseline_interval[1] & 
                                  as.vector(t(x[ , eval(time_column)])) <= baseline_interval[2]
                              }
                              # run through ids...
                              for (j in columns_used) {
                                if (use_data_table) {
                                  set(x = x, j = j, value = x[ , ..j] - as.numeric(baseline_means[ , .SD, .SDcols = j]) )
                                } else {
                                  baseline_mean <- mean(x[baseline_time_here, eval(j)], na.rm = TRUE)
                                  if (!is.na(baseline_mean)) {
                                    x[ , eval(j)] <- x[ , eval(j)] - baseline_mean
                                  } else {
                                    x[ , eval(j)] <- x[ , eval(j)] - mean(x[ , eval(j)], na.rm = TRUE)
                                  }
                                }
                              }
                              # return the updated subset
                              return(x)
                            }, 
                            columns_used = columns_used, time_column = time_column, baseline_interval = baseline_interval, 
                            use_data_table = use_data_table )
  
  # transform back to data table
  EEG_subset_new <- rbindlist(EEG_subset_list)
  rm(EEG_subset_list)
  
  # return the new edited copy
  return(EEG_subset_new)
}