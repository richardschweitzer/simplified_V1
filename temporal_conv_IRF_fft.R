# function to compute an fft-based IRF filter on the temporal dimension of a [y, x, t] matrix
# by Richard Schweitzer

# to debug:
# IRF = test_IRF$irf
# IRF_time = test_IRF$irf_time

temporal_conv_IRF_fft <- function(test_mat, IRF, IRF_time, 
                                  pad=1001, temporal_resolution=1000/1440, debug_mode=FALSE) {
  require(assertthat)
  require(torch)
  require(pracma)
  # a few assertions
  assert_that(length(dim(test_mat))==3, msg = "test_mat must have [y, x, t] dimension")
  assert_that(is.null(dim(IRF)), msg = "IRF must have one dimension")
  assert_that(mod(pad, 2)==1, msg = "pad must be odd!")
  assert_that(round(median(diff(IRF_time)), 4)==round(temporal_resolution, 4), 
              msg = "IRF does not have the same temporal resolution as specified by temporal_resolution!")
  # pad in a way that we get an even length of the sequence (allows not touching the DC component later)
  initial_len_seq <- test_mat$size(3)
  if (mod(initial_len_seq, 2)==1) { # sequence has odd length
    fft_len_seq <- initial_len_seq + pad 
  } else { # sequence has even length, let's not change that (pad is odd anyway)
    fft_len_seq <- initial_len_seq + pad + 1
  }
  # perform FFT here
  test_fft <- torch_fft_fft(self = test_mat, n = fft_len_seq, dim = 3) # 3 is temporal dimension
  if (debug_mode) {
    print(paste("dim(test_mat) =", paste(dim(test_mat), collapse = ",")))
    print(paste("dim(test_fft) =", paste(dim(test_fft), collapse = ",")))
  }
  # length of sequence
  len_seq <- test_fft$size(3)
  are_equal(fft_len_seq, len_seq)
  Fs <- round(1000/temporal_resolution)
  time_vec <- (0:(len_seq-1))*temporal_resolution   
  # get the frequencies and un-shift
  Fx <- Fs/len_seq*((-len_seq/2):(len_seq/2-1))
  Fx_shifted <- fftshift(Fx)
  # does the peak frequency match the intended frequency?
  if (debug_mode) {
    print(paste("Peak frequency is:", 
                paste(Fx_shifted[which.max(abs(as_array(test_fft[round(nrow(test_fft)/2), round(ncol(test_fft)/2), ])))], collapse = ",")
    ))
    plot(Fx_shifted, abs(as_array(test_fft[round(nrow(test_fft)/2), round(ncol(test_fft)/2), ])), type = "l")
  }
  # fft-transform the IRF
  IRF <- IRF / max(IRF)
  IRF_fft <- torch_fft_fft(self = IRF, n = fft_len_seq, dim = 1)
  # multiply FFT'ed matrix with FFT'ed IRF
  test_fft_irf <- torch_multiply(test_fft, IRF_fft) # this will be broadcasted
  are_equal(dim(test_fft_irf), dim(test_fft))
  if (debug_mode) {
    plot(Fx_shifted, abs(as_array(test_fft_irf[round(nrow(test_fft)/2), round(ncol(test_fft)/2), ])), 
         type = "l", col = "red")
    lines(Fx_shifted, abs(as_array(test_fft[round(nrow(test_fft)/2), round(ncol(test_fft)/2), ])))
  } else {
    rm(test_fft)
  }
  # inverse FFT transform
  test_mat_irf <- torch_fft_ifft(self = test_fft_irf, dim = 3)[ , , 1:initial_len_seq] # 3 is temporal dimension
  test_mat_irf <- test_mat_irf$real
  are_equal(dim(test_mat_irf), dim(test_mat))
  # do the sequences match?
  if (debug_mode) {
    plot(as_array(test_mat[round(nrow(test_mat)/2), round(ncol(test_mat)/2), ]))
    plot(as_array(test_mat_irf[round(nrow(test_mat)/2), round(ncol(test_mat)/2), ]))
  } else {
    rm(test_fft_irf, test_mat)
  }
  return(test_mat_irf)
}