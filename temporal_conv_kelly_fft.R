# function to compute an fft-based filter, according to Kelly's (1979) function, on the temporal dimension of a [y, x, t] matrix
# by Richard Schweitzer

temporal_conv_kelly_fft <- function(test_mat, SF, pad=1001, temporal_resolution=1000/1440, debug_mode=FALSE) {
  require(assertthat)
  require(torch)
  require(pracma)
  # a few assertions
  assert_that(length(dim(test_mat))==3)
  assert_that(exists("kelly_vel"))
  assert_that(mod(pad, 2)==1, msg = "pad must be odd!")
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
  # find Kelly sensitivity
  kelly_sens_Fx_shifted <- rep(0, length(Fx_shifted))
  kelly_sens_Fx_shifted[Fx_shifted!=0] <- kelly_vel(sf = SF, v = abs(Fx_shifted[Fx_shifted!=0]) / SF)
#  kelly_sens_Fx_shifted <- kelly_sens_Fx_shifted / max(kelly_sens_Fx_shifted)
  if (debug_mode) {
    plot(Fx_shifted, kelly_sens_Fx_shifted, type = "l")
  }
  are_equal(length(kelly_sens_Fx_shifted), len_seq)
  # multiply FFT'ed matrix with Kelly sensitivity
  test_fft_kelly <- torch_multiply(test_fft, kelly_sens_Fx_shifted) # this will be broadcasted
  are_equal(dim(test_fft_kelly), dim(test_fft))
  if (debug_mode) {
    plot(Fx_shifted, abs(as_array(test_fft_kelly[round(nrow(test_fft)/2), round(ncol(test_fft)/2), ])), 
         type = "l", col = "red")
    lines(Fx_shifted, abs(as_array(test_fft[round(nrow(test_fft)/2), round(ncol(test_fft)/2), ])))
  } else {
    rm(test_fft)
  }
  # inverse FFT transform
  test_mat_kelly <- torch_fft_ifft(self = test_fft_kelly, dim = 3)[ , , 1:initial_len_seq] # 3 is temporal dimension
  test_mat_kelly <- test_mat_kelly$real
  are_equal(dim(test_mat_kelly), dim(test_mat))
  # do the sequences match?
  if (debug_mode) {
    plot(as_array(test_mat[round(nrow(test_mat)/2), round(ncol(test_mat)/2), ]))
    plot(as_array(test_mat_kelly[round(nrow(test_mat)/2), round(ncol(test_mat)/2), ]))
  } else {
    rm(test_fft_kelly, test_mat)
  }
  return(test_mat_kelly)
}