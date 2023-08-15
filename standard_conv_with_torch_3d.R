standard_conv_with_torch_3d <- function(x, K, # x is a 3D torch tensor of size [1, t, y, x], K is a torch tensor of [y, x]
                                        no_CUDA=FALSE, use_half_on_GPU=FALSE ) {
  require(torch)
  require(assertthat)
  if (cuda_is_available()) {
    if (no_CUDA) {
      use_device <- torch_device('cpu')
    } else {
      use_device <- torch_device('cuda') # 'cuda'
    }
    # ... and what should be the preferred data type? (there are torch_half, torch_float, torch_double)
    if (use_half_on_GPU) {
      preferred_dtype <- torch_float16() # or half as mixed/half computing is optimized?
    } else {
      preferred_dtype <- torch_float32()
    }
  } else {
    use_device <- torch_device('cpu')
    preferred_dtype <- torch_float32()
  }
  # prepare the spatial kernel
  # note that pytorch implements a cross-correlation, not a convolution. 
  # Thus, we'll have opposite results as when done in the fft domain
  spat_kernel_1_gpu <- torch_tensor(K, device = use_device, dtype = preferred_dtype)
  spat_kernel_1_gpu <- torch_flip(spat_kernel_1_gpu, dims = c(1,2))
  # reshape spatial Gabor convolution kernel to (B, C, H, W)
  spat_kernel_1_gpu <- spat_kernel_1_gpu$reshape(c(1, 1, spat_kernel_1_gpu$size(1), spat_kernel_1_gpu$size(2)))
  spat_kernel_1_gpu <- spat_kernel_1_gpu$expand(c(dim(x)[2], 1, spat_kernel_1_gpu$size(3), spat_kernel_1_gpu$size(4)))
  # prepare the matrix, if it's already a CUDA GPU matrix, then don't do anything
  if (!check_if_torch_cuda(x)) {
    x_gpu <- torch_tensor(x, device = use_device, dtype = preferred_dtype)
  } else {
    x_gpu <- x
  }
  # run the spatial convolution, according to https://discuss.pytorch.org/t/manual-2d-convolution-per-channel/83907
  spatial_response_1 <- 
    torch_conv2d(input = x_gpu, weight = spat_kernel_1_gpu, 
                 # make sure the padding retains the overall size of the image
                 padding = c(spat_kernel_1_gpu$size(3)%/%2, spat_kernel_1_gpu$size(4)%/%2), 
                 groups = x_gpu$size(2) # number of time points
    )
  # return
  rm(spat_kernel_1_gpu, x_gpu)
  return(spatial_response_1)
}


### aux func:
check_if_torch_cuda <- function(x1) {
  # check's whether x is a torch CUDA object
  require(torch)
  x1_is_cuda <- FALSE
  try({ x1_is_cuda <- x1$is_cuda}, silent = TRUE)
  return(x1_is_cuda)
}
