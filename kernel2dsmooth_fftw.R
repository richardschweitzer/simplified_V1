kernel2dsmooth_fftw <- function (x, kernel.type = NULL, K = NULL, W = NULL, X = NULL, 
          xdim = NULL, Nxy = NULL, setup = FALSE, verbose = FALSE) 
  # modified for torch from the package 'smoothie'
{
  require(torch)
  if (is.null(xdim)) 
    xdim <- dim(x)
  if (is.null(Nxy)) 
    Nxy <- prod(xdim)
  if (is.null(W)) {
    if (!is.null(kernel.type)) 
      K <- kernel2dmeitsjer(type = kernel.type, ...)
    else if (is.null(K)) 
      stop("kernel2dsmooth: must give a value for at least one of kernel.type, K, or W")
    kdim <- dim(K)
    bigdim <- xdim + kdim - 1
    if (bigdim[1] <= 1024) 
      bigdim[1] <- 2^ceiling(log2(bigdim[1]))
    else bigdim[1] <- ceiling(bigdim[1]/512) * 512
    if (bigdim[2] <= 1024) 
      bigdim[2] <- 2^ceiling(log2(bigdim[2]))
    else bigdim[2] <- ceiling(bigdim[2]/512) * 512
    Kbig <- matrix(0, bigdim[1], bigdim[2])
    kcen <- floor((kdim + 1)/2)
    Kbig[1:(kdim[1] - kcen[1] + 1), 1:(kdim[2] - kcen[2] + 
                                         1)] <- K[kcen[1]:kdim[1], kcen[2]:kdim[2]]
    if (kdim[1] > 1) 
      Kbig[(bigdim[1] - kcen[1] + 2):bigdim[1], 1:(kdim[2] - 
                                                     kcen[2] + 1)] <- K[1:(kcen[1] - 1), kcen[2]:kdim[2]]
    if (kdim[2] > 1) 
      Kbig[1:(kdim[1] - kcen[1] + 1), (bigdim[2] - kcen[2] + 
                                         2):bigdim[2]] <- K[kcen[1]:kdim[1], 1:(kcen[2] - 
                                                                                  1)]
    if (all(kdim > 1)) 
      Kbig[(bigdim[1] - kcen[1] + 2):bigdim[1], (bigdim[2] - 
                                                   kcen[2] + 2):bigdim[2]] <- K[1:(kcen[1] - 1), 
                                                                                1:(kcen[2] - 1)]
    if (verbose) 
      cat("Finding the FFT of the kernel matrix.\n")
    W <- as_array(torch_fft_fft(torch_fft_fft(self = Kbig, 
                                              dim = 1), 
                                dim = 2)$cpu() )
    if (verbose) 
      cat("FFT of kernel matrix found.\n")
    if (setup) 
      return(W)
  }
  else bigdim <- dim(W)
  out <- matrix(0, bigdim[1], bigdim[2])
  out[1:xdim[1], 1:xdim[2]] <- x
  out[is.na(out)] <- 0
  if (verbose) 
    cat("Performing the convolution.\n")
  if (!is.null(X)) { 
    #out <- Re(fft(X * W, inverse = TRUE))[1:xdim[1], 1:xdim[2]]
    out <- torch_fft_ifft(torch_fft_ifft(self = torch_multiply(X, W), 
                                         dim = 1), 
                          dim = 2)[1:xdim[1], 1:xdim[2]]
    out <- as_array(out$real)
  } else { 
    #out <- Re(fft(fft(out) * W, inverse = TRUE))[1:xdim[1], 1:xdim[2]]
    out_fft <- torch_fft_fft(torch_fft_fft(self = out, 
                                           dim = 1), 
                             dim = 2)
    out_fft <- torch_multiply(out_fft, W)
    out <- torch_fft_ifft(torch_fft_ifft(self = out_fft, 
                                         dim = 1), 
                          dim = 2)[1:xdim[1], 1:xdim[2]]
    out <- as_array(out$real)
    rm(out_fft)
  }
  if (verbose) 
    cat("The convolution has been carried out.\n")
  return(zapsmall(out))
}