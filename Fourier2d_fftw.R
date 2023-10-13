Fourier2d_fftw <- function(x, bigdim = NULL, kdim = NULL) 
  # modified for torch from the package 'smoothie'
{
  require(torch)
  xdim <- dim(x)
  if (is.null(bigdim)) {
    if (is.null(kdim)) 
      stop("Fourier2d: one of bigdim or kdim must be supplied.")
    bigdim <- xdim + kdim - 1
    if (bigdim[1] <= 1024) 
      bigdim[1] <- 2^ceiling(log2(bigdim[1]))
    else bigdim[1] <- ceiling(bigdim[1]/512) * 512
    if (bigdim[2] <= 1024) 
      bigdim[2] <- 2^ceiling(log2(bigdim[2]))
    else bigdim[2] <- ceiling(bigdim[2]/512) * 512
  }
  out <- matrix(0, bigdim[1], bigdim[2])
  out[1:xdim[1], 1:xdim[2]] <- x
  out[is.na(out)] <- 0
  fft_out <- as_array(torch_fft_fft(torch_fft_fft(self = out, 
                                                  dim = 1), 
                                    dim = 2)$cpu() )
  return(fft_out)
}
