fft_conv_with_torch_3d <- function(x, K, # x is a 3D torch tensor of size [1, t, y, x], K is a torch tensor of [y, x]
                                   use_torch=FALSE, no_CUDA=FALSE, use_half_precision=FALSE, 
                                   debug_mode=FALSE, 
                                   use_compact_writing=TRUE # to see whether this makes a difference in speed
                                   ) {
  # based on kernel2dsmooth from package 'smoothie', changed to iterate over time slices
  require(torch)
  require(assertthat)
  if (use_torch) {
    if (cuda_is_available()) {
      if (no_CUDA) {
        use_device <- torch_device('cpu')
      } else {
        use_device <- torch_device('cuda') # 'cuda'
      }
    } else {
      use_device <- torch_device('cpu')
    }
    # ... and what should be the preferred data type? (there are torch_half, torch_float, torch_double)
    if (use_half_precision) {
      preferred_dtype <- torch_float16() 
    } else {
      preferred_dtype <- torch_float32()
    }
  }
  # IMPORTANT: torch_fft only works with torch_float32 !!!
  
  # get dimensions for padding right
  xdim <- dim(x)[3:4] # dimension of image, ignoring batch and channel, e.g., 1 482 167 254
  assert_that(dim(x)[1]==1, msg = "Batch number must be one!")
  Nxy <- prod(xdim)
  kdim <- dim(K)
  bigdim <- xdim + kdim - 1
  if (debug_mode) {
    print(paste("xdim:", paste(xdim, collapse = ",")))
    print(paste("kdim:", paste(kdim, collapse = ",")))
  }
  if (bigdim[1] <= 1024) {
    bigdim[1] <- 2^ceiling(log2(bigdim[1]))
  } else {
    bigdim[1] <- ceiling(bigdim[1]/512) * 512
  }
  if (bigdim[2] <= 1024) {
    bigdim[2] <- 2^ceiling(log2(bigdim[2]))
  } else {
    bigdim[2] <- ceiling(bigdim[2]/512) * 512
  }
  if (debug_mode) {
    print(paste("bigdim:", paste(bigdim, collapse = ",")))
  }
  # create padded Kernel matrix
  if (use_torch) {
    Kbig <- torch_zeros(c(bigdim[1], bigdim[2]), 
                        device = use_device, dtype = torch_float32()) # must be full precision for fft
  } else {
    Kbig <- matrix(0, bigdim[1], bigdim[2])
  }
  kcen <- floor((kdim + 1)/2)
  Kbig[1:(kdim[1] - kcen[1] + 1), 
       1:(kdim[2] - kcen[2] + 1)] <- K[kcen[1]:kdim[1], 
                                       kcen[2]:kdim[2]]
  if (kdim[1] > 1) {
    Kbig[(bigdim[1] - kcen[1] + 2):bigdim[1], 
         1:(kdim[2] - kcen[2] + 1)] <- K[1:(kcen[1] - 1), 
                                         kcen[2]:kdim[2]]
  }
  if (kdim[2] > 1) {
    Kbig[1:(kdim[1] - kcen[1] + 1), 
         (bigdim[2] - kcen[2] + 2):bigdim[2]] <- K[kcen[1]:kdim[1], 
                                                   1:(kcen[2] - 1)]
  }
  if (all(kdim > 1)) {
    Kbig[(bigdim[1] - kcen[1] + 2):bigdim[1], 
         (bigdim[2] - kcen[2] + 2):bigdim[2]] <- K[1:(kcen[1] - 1), 
                                                   1:(kcen[2] - 1)]
  }
  # perform FFT of kernel matrix
  if (use_torch) {
    W <- torch_fft_fft(torch_fft_fft(self = Kbig, 
                                     dim = 1), 
                       dim = 2) #/ prod(bigdim) # this division is not needed in torch
  } else {
    W <- fft(Kbig)/prod(bigdim)
  }
  if (debug_mode) {
    print(paste("dim(W):", paste(dim(W), collapse = ",")))
  }
  rm(Kbig, K) # don't need it anymore, as we now have it's Fourier transform
  # prepare the target matrix
  if (use_torch) {
    out <- torch_zeros(c(dim(x)[1], dim(x)[2], bigdim[1], bigdim[2]), 
                       device = use_device, dtype = torch_float32()) # must be full precision for fft
  } else {
    out <- array(data = 0, dim = c(dim(x)[1], dim(x)[2], bigdim[1], bigdim[2]))
  }
  if (debug_mode) {
    print(paste("padded dim(out):", paste(dim(out), collapse = ",")))
  }
  out[1, 1:(dim(x)[2]), 1:(xdim[1]), 1:(xdim[2])] <- x
  rm(x) # is not needed anymore as now part of 'out'
  if (use_torch) {
    out[torch_isnan(out)==1] <- 0
  } else {
    out[is.na(out)] <- 0
  }
  # now walk through time slices
  for (t_i in 1:(dim(out)[2])) {
    # if (debug_mode) {
    #   print(t_i)
    # }
    if (use_torch) {
      if (torch_is_nonzero(torch_sum(torch_abs(out[1, t_i, , ])))) { # check for nonzeros
        if (!use_compact_writing) {
          out_fft <- torch_fft_fft(torch_fft_fft(self = out[1, t_i, , ], 
                                                 dim = 1), 
                                   dim = 2)
          out_fft <- torch_multiply(out_fft, W)
          out_ifft <- torch_fft_ifft(torch_fft_ifft(self = out_fft, 
                                                    dim = 1), 
                                     dim = 2)
        } else {
          out_ifft <- 
            torch_fft_ifft(
              torch_fft_ifft(
                torch_multiply(
                  torch_fft_fft(
                    torch_fft_fft(
                      out[1, t_i, , ], 
                      dim = 1), 
                    dim = 2), 
                  W), 
                dim = 1), 
              dim = 2)
        }
        out[1, t_i, , ] <- out_ifft$real
      } 
    } else {
      if (any(out[1, t_i, , ] != 0)) { # check for nonzeros
        out[1, t_i, , ] <- Re(fft(fft(out[1, t_i, , ]) * W, inverse = TRUE))
      }
    }
  }
  # clean up
  rm(W)
  if (use_torch) { 
    if (use_compact_writing) {
      rm(out_ifft)
    } else {
      rm(out_fft, out_ifft) 
    }
  }
  # and reduce to original size
  out <- out[1:(dim(out)[1]), 1:(dim(out)[2]), 1:xdim[1], 1:xdim[2]]
  if (debug_mode) {
    print(paste("unpadded dim(out):", paste(dim(out), collapse = ",")))
  }
  # if necessary, roll back to half precision
  if (preferred_dtype == torch_float16()) {
    out <- out$to(dtype = preferred_dtype, device = use_device)
  }
  return(out)
}


