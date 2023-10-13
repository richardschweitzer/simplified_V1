pad_and_resample <- function(x, p, q, d=10, perc = 2/10) {
  require(signal)
  # do we downsample? --> will resample's FIR filter be used?
  if (p/q < 1) {
    # compute the delay to compensate for FIR filter
    # see: https://de.mathworks.com/help/signal/ug/compensate-for-the-delay-introduced-by-an-fir-filter.html
    # and: https://dsp.stackexchange.com/questions/18435/group-delay-of-the-fir-filter
    # N, the filter length, is n=2*d+1 in resample, thus the group delay of the FIR is (n-1)/2, i.e., d
    delay = ceiling((d*2+1)/2) # from Matlab's resample
  } else {
    delay = 0
  }
  # how many pads?
  n_pads <- ceiling(length(x) * perc)
  # median of the signal to normalize
  m_x <- median(x, na.rm = TRUE)
  # apply the padding, apply the predicted delay
  x_padded <- c(rep(x = x[1], times = n_pads-delay), 
                x, 
                rep(x = x[length(x)], times = n_pads+delay)) - m_x
  # run the resampling
  future_size_of_x <- round(length(x) * (p / q))
  x_new_padded <- resample(x = x_padded, p = p, q = q, d = d)
  # how many of the samples of the padding to remove?
  n_pads_new <- round(n_pads * (length(x_new_padded)/length(x_padded)))
  # remove the padding
  x_new <- x_new_padded[(n_pads_new+1):(n_pads_new+future_size_of_x)] + m_x
#  x_new <- x_new_padded[(n_pads_new+1):(length(x_new_padded)-n_pads_new)] + m_x
  return(x_new)
}


resample_signal_df <- function(df, time_col, id_cols, linear_cols, signal_cols, constant_cols, new_samp_rate, 
                               plot_results=FALSE, debug_mode=FALSE) {
  # This function resamples a data.frame df of e.g. eye movement data
  # - time_col: column where timestamps are placed [in milliseconds]
  # - id_cols: columns where no interpolation should take place, as values are constant or non-numeric IDs
  # - linear_cols: columns where linear interpolation should take place, as values are metric
  # - spline_cols: columns where spline interpolation should take place, as values are metric
  # - constant_cols: columns where constant interpolation should take place, e.g., when values are 0 and 1
  # - new_samp_rate: sampling rate [Hz] we want to resample to. 
  # - plot_results: optional plot of all the interpolated values (default: FALSE)
  # by Richard Schweitzer, 07/2019
  
  # debug with:
  # df = eye_sac # (from STRK analyses)
  # time_col = "sac_t_raw"
  # id_cols = c("session_name", "subj_id", "session_id", "trial_id", "ID", "stim_present", "nr_iteration", "experiment")
  # linear_cols = c("sac_t_pc", "sac_x_raw", "sac_y_raw", "sac_pupil")
  # signal_cols = c("sac_t_pc", "sac_x_raw", "sac_y_raw", "sac_pupil") # , "sac_x_raw_2", "sac_y_raw_2", "sac_pupil_2"
  # constant_cols = c()
  # new_samp_rate = 250
  # plot_results = TRUE
  
  require(assertthat)
  require(data.table)
  
  # check data.frame
  df <- as.data.frame(df)
  assert_that(nrow(df)>5) # we don't need data we cannot interpolate from
  
  # check whether id_cols are all unique
  unique_id_cols <- sapply(df[ , is.element(colnames(df), id_cols)], FUN = unique)
  if (debug_mode) print(unique_id_cols)
  assert_that(all(sapply(unique_id_cols, length)==1)) # id cols must contain unique values 
  
  # check if time series is sorted
  timestamps <- df[ , is.element(colnames(df), time_col)]
  assert_that(!is.unsorted(timestamps))
  
  # compute p and q based on the current and desired sampling rate
  current_samp_rate <- 1000/median(diff(timestamps), na.rm = TRUE)
  if (current_samp_rate >= new_samp_rate) {
    p = 1
    q = current_samp_rate / new_samp_rate
  } else {
    p = new_samp_rate / current_samp_rate
    q = 1
  }
  target_size = round(length(timestamps) * (p / q)) # or 'ceiling', as described in the help
  
  # make new timestamps according to the new_samp_rate
  new_timestamps <- seq(from = min(timestamps), to = max(timestamps), by = 1000/new_samp_rate)[1:target_size]
  
  # here we create the new df
  new_df <- data.frame(new_timestamps)
  colnames(new_df) <- time_col
  
  # LINEAR interpolation here
  for (co in linear_cols) { # for all linear_cols, do this:
    original_data <- df[ , is.element(colnames(df), co)]
    if (length(original_data[!is.na(original_data)])>2) {
      func <- approxfun(x = timestamps, y = original_data, 
                        method = "linear") # fit the linear interpolation function
      temp_data <- func(new_timestamps)
    } else {
      temp_data <- rep(NaN, length(new_timestamps))
    }
    new_data <- data.frame(temp_data) # interpolate
    colnames(new_data) <- co
    new_df <- cbind(new_df, new_data)
  }
  
  # NEW: bandlimited interpolation here
  for (co in signal_cols) { # for all spline_cols, do this:
    original_data <- df[ , is.element(colnames(df), co)]
    resampled_data <- pad_and_resample(x = original_data, p = p, q = q) # a <- resample(x = original_data, p = p, q = q, d = 10)
    new_data_resampled <- data.frame(resampled_data) # interpolate
    colnames(new_data_resampled) <- co
    new_df <- cbind(new_df, new_data_resampled)
  }
  # if (debug_mode) {
  #   plot(timestamps[1:100], original_data[1:100]) #, ylim = c(min(temp_data), max(temp_data))
  #   lines(new_timestamps[1:50], temp_data[1:50], col = "blue")
  #   lines(new_timestamps[1:50], resampled_data[1:50], col = "red")
  #   legend("bottomleft", c("Original", "approxfun (base)", "resample (signal)"),
  #          col = c("black", "blue", "red"),
  #          pch = c(1, NA, NA),
  #          lty = c(NA, 1, 1), bty = "n")
  # }
  
  # CONSTANT interpolation here
  for (co in constant_cols) { # for all constant_cols, do this:
    original_data <- df[ , is.element(colnames(df), co)]
    func <- approxfun(x = timestamps, y = original_data, 
                      method = "constant") # fit the linear interpolation function
    new_data <- data.frame(func(new_timestamps)) # interpolate
    colnames(new_data) <- co
    new_df <- cbind(new_df, new_data)
  }
  
  # ID cols are added now... (we just use the first row, as they should be all the same)
  id_df <- df[1, is.element(colnames(df), id_cols)]
  id_df <- rbindlist(replicate(nrow(new_df), id_df, simplify = FALSE))
  new_df <- cbind(new_df, id_df) 
  rm(id_df)
  
  # plot results, if this is desired:
  if (plot_results) {
    # first, the linear cols
    if (length(linear_cols)>0) {
      par(mfrow = c(ceiling((length(linear_cols))/2), 2))
      par(mar=c(4,4,4,4))
      for (co in linear_cols) {
        plot(df[ ,is.element(colnames(df), time_col)], 
             df[ ,is.element(colnames(df), co)], 
             xlab = "time", ylab = co)
        points(new_df[ ,is.element(colnames(new_df), time_col)], 
               new_df[ ,is.element(colnames(new_df), co)], col = 2, pch = "*")
      }
    }
    # second, the resample cols
    if (length(signal_cols)>0) {
      par(mfrow = c(ceiling((length(signal_cols))/2), 2))
      par(mar=c(4,4,4,4))
      for (co in signal_cols) {
        plot(df[ ,is.element(colnames(df), time_col)], 
             df[ ,is.element(colnames(df), co)], 
             xlab = "time", ylab = co)
        points(new_df[ ,is.element(colnames(new_df), time_col)], 
               new_df[ ,is.element(colnames(new_df), co)], col = 2, pch = "*")
      }
    }
    # third, the constant cols
    if (length(constant_cols)>0) {
      par(mfrow = c(ceiling((length(constant_cols))/2), 2))
      par(mar=c(4,4,4,4))
      for (co in constant_cols) {
        plot(df[ ,is.element(colnames(df), time_col)], 
             df[ ,is.element(colnames(df), co)], 
             xlab = "time", ylab = co)
        points(new_df[ ,is.element(colnames(new_df), time_col)], 
               new_df[ ,is.element(colnames(new_df), co)], col = 2, pch = "*")
      }
    }
    # reset
    par(mfrow = c(1,1))
  }
  
  # final product:
  return(new_df)
}


# # get someexample data here
# df = trial_eye
# new_samp_rate = 500
# time_col = "time_DP"
# id_cols = c("ID", "experiment", "trial_id", "session_id", "subj_id", "nr_iteration")
# linear_cols = c("time_pc", "right_x", "right_y", "left_x", "left_y", "right_pupil", "left_pupil" )
# constant_cols = c("right_blink", "left_blink")
# 
# # example run: 
# new_df <- resample_df(df, time_col, id_cols, linear_cols, constant_cols, 500, TRUE)
