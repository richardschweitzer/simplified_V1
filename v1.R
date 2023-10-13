v1 <- function(mat_over_t = NULL, # you can pass a full stimulus over time matrix [y,x,t] to the function
               stim_mat, # the noise patch as a matrix
               gabor_list, # a list of gabor filters created by 'get_gabor_filter_bank'
               normalize_RFs = TRUE, # TRUE: to sum-divide the Gabor filter kernels
               irf_df = NULL, # a data.frame of IRFs resolved in SF x time created by 'get_IRF_df'
               normalize_IRFs = FALSE, # TRUE: to sum-divide the IRFs
               signal_x, signal_y, signal_t, # the signal in pixels
               norm_sigma = NULL, # provide the normalization constant, if not it will be estimated based on the max
               use_normalization_pool = TRUE, # TRUE: use a normalization pool in space, SF, and Orientation
               normalize_override = FALSE, # TRUE: to simply not do normalization
               output_full_sequences = FALSE, # TRUE: returns complete sequences. WARNING: this is memory-intense
               output_full_sequences_use_RleArray = FALSE, # TRUE: use RleArray (from package DelayedArray) if output_full_sequences=TRUE
               signal_x_range = NULL, signal_y_range = NULL, signal_t_range = NULL, # the min and max of temporal scale
               no_CUDA = FALSE, # set to TRUE if you want to use CPU instead of GPU (but torch)
               use_half_precision = FALSE, # if a GPU is available, use half-precision?
               spatial_resolution = 1, # spatial resolution of processing function, i.e., in pixels
               temporal_resolution = 1000/1440, # temporal resolution of processing function in milliseconds
               debug_mode = FALSE, # shows output at every step
               show_final_maps = TRUE,  # shows all resulting 2D maps at the end
               final_maps_filename = "v1_output.pdf" # name of the pdf file if show_final_maps==TRUE
) {
  options(torch.threshold_call_gc = 4000)
  require(torch)
  require(data.table)
  require(assertthat)
  require(pracma)
  require(tictoc)
  if (output_full_sequences & output_full_sequences_use_RleArray) {
    # to deal with massive objects in memory, required: install.packages("BiocManager")
    # then: BiocManager::install("DelayedArray")
    require(DelayedArray)
  }
  ## 0.1 quick checks
  # do the SFs of the supplied IRF and Gabor objects match?
  gabor_SFs <- sort(unique(as.numeric(unlist(lapply(X = gabor_list, FUN = function(x) { x[[3]]['rf_freq_dva'] })))))
  if (!is.null(irf_df)) {
    irf_SFs <- sort(unique(irf_df$SF))
    are_equal(gabor_SFs, irf_SFs)
  }
  # have we loaded the convolution functions?
  assert_that(exists("standard_conv_with_torch_3d"))
  assert_that(exists("fft_conv_with_torch_3d"))
  assert_that(exists("plot_heatmap"))
  if (is.null(irf_df)) { # no IRFs provided, thus use Kelly-filter approach
    assert_that(exists("kelly_vel"))
    assert_that(exists("temporal_conv_kelly_fft"))
  }
  # is the IRF properly scaled in the temporal domain?
  if (!is.null(irf_df)) { # IRFs provided
    irf_Time <- sort(unique(irf_df$irf_time))
    are_equal(length(unique(round(diff(irf_Time), 5))), 1)
    are_equal(unique(round(diff(irf_Time), 5)), round(temporal_resolution, 5))
  }
  ## 0.2 compute properties of signal
  if (is.null(mat_over_t)) { # construct stimulus matrix
    length_signal_x <- abs(max(signal_x)-min(signal_x)) # the spatial extent of the signal 
    length_signal_y <- abs(max(signal_y)-min(signal_y))
    length_signal_t <- abs(max(signal_t)-min(signal_t)) # the spatial extent of the signal 
    assert_that(length(signal_x)==length(signal_y))
    assert_that(length(signal_x)==length(signal_t))
    n_samples <- length(signal_t)
    signal_df <- data.frame(signal_x, signal_y, signal_t)
    signal_df$n_samples <- 1:n_samples
    if (debug_mode) {
      print("Input:")
      print(signal_df)
    }
  } else { # stimulus matrix is already there
    if (debug_mode) {
      print("You have provided a mat_over_t. Its dimensions:")
      print(dim(mat_over_t))
    }
  }
  ## STEP 1: fill the position-over-time matrix
  if (is.null(mat_over_t)) { # indeed, we have to construct the STIMULUS MATRIX
    
    # 1.1 get the dimensions for the cells of the matrix
    spatial_pad_size <- round(max(dim(stim_mat))) # times two ???
    if (!is.null(irf_df)) { # IRFs provided
      temporal_pad_size <- round(1*(max(irf_Time)-min(irf_Time)))
    } else {
      temporal_pad_size <- 300
    }
    # X
    if (!is.null(signal_x_range)) {
      mat_size_x = seq(signal_x_range[1]-spatial_pad_size, signal_x_range[2]+spatial_pad_size, 
                       by = spatial_resolution) # columns, x coordinate
    } else {
      mat_size_x = seq(min(signal_x)-spatial_pad_size, max(signal_x)+spatial_pad_size, 
                       by = spatial_resolution) # columns, x coordinate
    }
    # Y
    if (!is.null(signal_y_range)) {
      mat_size_y = seq(signal_y_range[1]-spatial_pad_size, signal_y_range[2]+spatial_pad_size, 
                       by = spatial_resolution) # rows, y coordinate
    } else {
      mat_size_y = seq(min(signal_y)-spatial_pad_size, max(signal_y)+spatial_pad_size, 
                       by = spatial_resolution) # rows, y coordinate
    }
    # T
    if (!is.null(signal_t_range)) { # compute the temporal scale based on signal_t_range
      mat_size_t = seq(signal_t_range[1], signal_t_range[2]+temporal_pad_size, 
                       by = temporal_resolution)
    } else { # compute the temporal scale based on signal_t
      min_t <- min(signal_t) # should this be fixed at zero?
      mat_size_t = seq(min_t, max(signal_t)+temporal_pad_size, 
                       by = temporal_resolution) # third dimension, t coordinate (0 was previously: min(signal_t)-round(length_signal_t))
    }
    # 1.3 create large matrix of zeros according to the specified resolution
    mat_over_t <- array(data = 0, dim = c(length(mat_size_y), length(mat_size_x), length(mat_size_t)), 
                        dimnames = list(mat_size_y, mat_size_x, mat_size_t) )
    print(paste("dim(mat_over_t) =", paste0(dim(mat_over_t), collapse = ","), "[y,x,t]" ))
    # 1.4 fill the matrix with the blob where the position is present in the signal
    if (debug_mode) {
      p_stim_sequence_list <- vector(mode = "list", length = n_samples)
    }
    for (sample_i in (1:n_samples)) {
      # in case we have an NA, we'll leave it out
      if (!is.na(signal_x[sample_i]) && !is.na(signal_y[sample_i]) && !is.na(signal_t[sample_i])) {
        put_sample_x_here <- which.min(abs(mat_size_x - signal_x[sample_i]))
        put_sample_y_here <- which.min(abs(mat_size_y - signal_y[sample_i]))
        put_sample_t_here <- which.min(abs(mat_size_t - signal_t[sample_i]))
        mat_over_t[seq(put_sample_y_here-round(nrow(stim_mat)/2), 
                       put_sample_y_here-round(nrow(stim_mat)/2)+nrow(stim_mat)-1, 
                       by = 1), 
                   seq(put_sample_x_here-round(ncol(stim_mat)/2), 
                       put_sample_x_here-round(ncol(stim_mat)/2)+ncol(stim_mat)-1, 
                       by = 1), 
                   put_sample_t_here] <- stim_mat
        #print(address(mat_over_t)) # this operation works in place
        if (debug_mode) {
          p_stim_now <- plot_heatmap(mat_over_t[ , , put_sample_t_here], do_print = FALSE)
          p_stim_sequence_list[[sample_i]] <- p_stim_now + 
            #ggtitle(paste("t =", signal_t[sample_i])) + 
            theme(legend.position = "none")
        }
      }
    }
    # plot the sequence?
    if (debug_mode) {
      require(cowplot)
      p_stim_sequence <- plot_grid(plotlist = p_stim_sequence_list)
      print(p_stim_sequence)
    }
    last_signal_t_index <- put_sample_t_here # save the last indeces where we have a stimulus present
    
  } else { # the matrix is already provided
    
    if (is.null(irf_df)) {
      temporal_pad_size <- 300
    }
    mat_size_y = as.numeric(dimnames(mat_over_t)[[1]])
    mat_size_x = as.numeric(dimnames(mat_over_t)[[2]])
    mat_size_t = as.numeric(dimnames(mat_over_t)[[3]])
    last_signal_t_index <- min(which(apply(X = mat_over_t, 
                                           FUN = function(x) {all(x==0)}, 
                                           MARGIN = c(3)))) - 1
    if (is.infinite(last_signal_t_index) || is.null(last_signal_t_index)) {
      last_signal_t_index <- round(length(mat_size_t)/2)
    }
    
  } # stimulus matrix done.
  
  ## STEP 2: Convolutions performed in loop over SF and Ori
  # is cuda available?
  if (cuda_is_available()) {
    if (debug_mode) {
      print(paste("CUDA is available! Version:", cuda_runtime_version()))
    }
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
    preferred_dtype <- torch_float16() # or half as mixed/half computing is optimized?
    if (use_device == torch_device('cpu')) { # sorry, on the CPU most operations are only allowed with full precision
      use_half_precision <- FALSE
      preferred_dtype <- torch_float32()
    }
  } else {
    preferred_dtype <- torch_float32()
  }
  # 0) reshape the mat_over_t to match the 4D requirements of the 2D convolution
  # transfer to GPU, if necessary:
  mat_over_t_gpu <- torch_tensor(mat_over_t, device = use_device, dtype = preferred_dtype)
  mat_over_t_gpu <- mat_over_t_gpu$permute(c(3, 1, 2))$unsqueeze(1)
  if (debug_mode) {
    print("Converted matrix to tensor and reshaped to dimensions:")
    print(mat_over_t_gpu$size())
  }
  # here we store SFs and Oris for each filter
  SFs_of_filters <- vector(mode = "numeric", length = length(gabor_list))
  Oris_of_filters <- vector(mode = "numeric", length = length(gabor_list))
  # here we store squared and low-pass responses for each filter
  lowpass_3d_over_filters <- vector(mode = "list", length = length(gabor_list))
  squared_3d_over_filters <- vector(mode = "list", length = length(gabor_list))
  # this is THE LOOP to compute squared and lowpass responses over filters:
  for (gabor_list_i in 1:length(gabor_list)) {
    # which SF are we dealing with?
    SF_now <- gabor_list[[gabor_list_i]][[3]]['rf_freq_dva']
    Ori_now <- gabor_list[[gabor_list_i]][[3]]['rf_ori']
    if (debug_mode || show_final_maps) { 
      print(paste0("Gabor filter ", gabor_list_i, " of ", length(gabor_list), 
                   ", SF=", round(SF_now, 2), ", Ori=", round(Ori_now, 2) )) 
    }
    # get corresponding IRF and normalize
    if (!is.null(irf_df)) { # IRFs provided
      IRF_now <- irf_df[SF==SF_now, irf]
      IRF_now_time <- irf_df[SF==SF_now, irf_time]
      if (normalize_IRFs) {
        IRF_now <- IRF_now / sum(abs(IRF_now))
      }
    } else {
      IRF_now_time <- seq(0, temporal_pad_size, by = temporal_resolution)
    }
    # get Iris Groen's low-pass filter for the delayed-normalization model
    prm.tau2 <- 0.75 # see Figure 11 from "Temporal Dynamics of Neural Responses in Human Visual Cortex"
    IRF_lp = exp(-(IRF_now_time/1000)/prm.tau2)
    if (normalize_IRFs) {
      IRF_lp = IRF_lp / sum(abs(IRF_lp)) 
    }
    # get Gabor filter kernels (1: zero, 2: half phase) and make sure they have odd dimensions
    spat_kernel_1 <- gabor_list[[gabor_list_i]][[1]]
    spat_kernel_2 <- gabor_list[[gabor_list_i]][[2]]
    assert_that(all(mod(dim(spat_kernel_1), 2) == 1), msg = "spat_kernel_1 does not have odd dimensions")
    assert_that(all(mod(dim(spat_kernel_2), 2) == 1), msg = "spat_kernel_2 does not have odd dimensions")
    # normalize
    if (normalize_RFs) {
      spat_kernel_1 <- spat_kernel_1 / sum(abs(spat_kernel_1))
      spat_kernel_2 <- spat_kernel_2 / sum(abs(spat_kernel_2))
    }
    ## 1) determine how convolution should be performed - and perform it
    # compute number of operations. If there are more than 10^8.5 in the spatial domain, use FFT-based
    number_of_spatial_operations <- prod(dim(spat_kernel_1)) * 
      mat_over_t_gpu$size(3) * mat_over_t_gpu$size(4)
    if ( (cuda_is_available() && number_of_spatial_operations > 10^8.5 #10^9 
          ) ||
         ((!cuda_is_available() || no_CUDA) && number_of_spatial_operations > 10^7
          )
    ) {
      use_fft_based_2dconv <- TRUE
    } else {
      use_fft_based_2dconv <- FALSE
    }
    ## 2) run the spatial convolution
    # note that pytorch implements a cross-correlation, not a convolution. 
    # Thus, we'll have opposite results as when done in the fft domain
    spatial_response_1 <- NULL
    spatial_response_2 <- NULL
    if (use_fft_based_2dconv) { # fft-based convolution (does not work with half precision!)
      tryCatch({ # to try to catch the case where we run out of CUDA memory
        tic(msg = "fft-based 2D convolution")
        no_CUDA_for_fft_conv <- FALSE # we may not have enough memory to do this on the GPU...
        spatial_response_1 <- fft_conv_with_torch_3d(x = mat_over_t_gpu, K = spat_kernel_1, use_torch = TRUE,
                                                     no_CUDA = no_CUDA_for_fft_conv, use_half_precision = use_half_precision)
        spatial_response_2 <- fft_conv_with_torch_3d(x = mat_over_t_gpu, K = spat_kernel_2, use_torch = TRUE,
                                                     no_CUDA = no_CUDA_for_fft_conv, use_half_precision = use_half_precision)
        # in case we used the CPU we might have to take those spatial responses back to the GPU
        if (no_CUDA_for_fft_conv==TRUE && no_CUDA==FALSE && cuda_is_available() ) {
          spatial_response_1 <- spatial_response_1$to(device = use_device, dtype = preferred_dtype)
          spatial_response_2 <- spatial_response_2$to(device = use_device, dtype = preferred_dtype)
        }
        toc(quiet = !(debug_mode || show_final_maps))
      }, error = function(e) { print(substr(e$message, 1, 250)) } )
    }
    if (is.null(spatial_response_1) | is.null(spatial_response_2)) { # normal convolution, also in case the fft-based one did not work
      tic(msg = "standard 2D convolution")
      spatial_response_1 <- standard_conv_with_torch_3d(x = mat_over_t_gpu, K = spat_kernel_1, 
                                                        no_CUDA = no_CUDA, use_half_precision = use_half_precision)
      spatial_response_2 <- standard_conv_with_torch_3d(x = mat_over_t_gpu, K = spat_kernel_2, 
                                                        no_CUDA = no_CUDA, use_half_precision = use_half_precision)
      toc(quiet = !(debug_mode || show_final_maps))
    }
    assert_that(all(mat_over_t_gpu$size()==spatial_response_1$size()))
    assert_that(all(mat_over_t_gpu$size()==spatial_response_2$size()))
    # ... diagnostic plots: what does the convolution look like?
    if (debug_mode) {
      print("Performed 2D spatial convolution, thereby achieving matrix with size:")
      print(paste(paste0(spatial_response_1$size(), collapse = ","), "and", paste0(spatial_response_2$size(), collapse = ",")))
      p_spatial_responses <- vector(mode = "list", length = 5)
      p_spatial_responses[[1]] <- plot_heatmap(as_array(mat_over_t_gpu[1, last_signal_t_index-3, , ]$cpu()), do_print = FALSE) + 
        ggtitle("Original stimulus")
      both_spatial_kernels <- rbind(spat_kernel_1, spat_kernel_2)
      dimnames(both_spatial_kernels) <- NULL
      p_spatial_responses[[2]] <- plot_heatmap(both_spatial_kernels, do_print = FALSE) + 
        ggtitle("Spatial kernels (zero & half)")
      p_spatial_responses[[3]] <- plot_heatmap(as_array(spatial_response_1[1, last_signal_t_index-3, , ]$cpu()), do_print = FALSE) + 
        ggtitle(paste0("Spatial resp (zero), SF=", round(SF_now, 2), ", Ori=", round(Ori_now, 2) ))
      p_spatial_responses[[4]] <- plot_heatmap(as_array(spatial_response_2[1, last_signal_t_index-3, , ]$cpu()), do_print = FALSE) + 
        ggtitle(paste0("Spatial resp (half), SF=", round(SF_now, 2), ", Ori=", round(Ori_now, 2) ))
      p_spatial_responses[[5]] <- plot_heatmap(sqrt(as_array(spatial_response_1[1, last_signal_t_index-3, , ]$cpu())^2 + 
                                                      as_array(spatial_response_2[1, last_signal_t_index-3, , ]$cpu())^2), do_print = FALSE) + 
        ggtitle(paste0("Squared responses"))
      p_spatial_responses <- plot_grid(plotlist = p_spatial_responses, nrow = 3)
      print(p_spatial_responses)
      rm(both_spatial_kernels)
    } else {
      rm(spat_kernel_1, spat_kernel_2)  # save memory!
    }
    ## 3) prepare the temporal response kernel (is has to have odd dimension - and the procedure below assures that)
    # Note that this will only work with full precision
    # spatial_response is now [1, t, y, x]
    if (!is.null(irf_df)) { # IRFs provided
      resp_fun_gpu <- torch_tensor(data = IRF_now, device = use_device, dtype = preferred_dtype)
      resp_fun_gpu <- resp_fun_gpu$expand(c(spatial_response_1$size(4), 1, resp_fun_gpu$size())) # expand to match matrix dimensions (x dimension)
      resp_fun_gpu <- torch_cat(list(torch_zeros(c(resp_fun_gpu$size(1), resp_fun_gpu$size(2), resp_fun_gpu$size(3)-1), 
                                                 device = use_device, dtype = preferred_dtype), # create the zero pad to center the TRF
                                     resp_fun_gpu), dim = 3)
      if (debug_mode) {
        print("Created temporal response tensor of dimensions for 1D convolution:")
        print(resp_fun_gpu$size())
      }
    } else {
      if (debug_mode) {
        print("No temporal response function provided, thus using FFT-based Kelly function.")
      }
    }
    ## 4) prepare and run the temporal convolution
    tic(msg = "Convolutions through time")
    spatial_response_1 <- spatial_response_1$permute(c(1, 3, 4, 2))$squeeze() # new dim: y, x, t, as initially
    spatial_response_2 <- spatial_response_2$permute(c(1, 3, 4, 2))$squeeze()
    if (debug_mode) {
      print("Reshaped main tensor to dimensions to prepare for 1D temporal convolution:")
      print(spatial_response_1$size())
    }
    if (!is.null(irf_df)) { # IRFs provided
      temporal_response_1 <- torch_conv1d(input = spatial_response_1, 
                                          weight = resp_fun_gpu$flip(3), 
                                          padding = resp_fun_gpu$size(3)%/%2, 
                                          groups = spatial_response_1$size(2) # is the x dimension
      )
      temporal_response_2 <- torch_conv1d(input = spatial_response_2, 
                                          weight = resp_fun_gpu$flip(3), 
                                          padding = resp_fun_gpu$size(3)%/%2, 
                                          groups = spatial_response_2$size(2)
      )
    } else { # use Kelly's function directly in an FFT-based approach
      temporal_response_1 <- temporal_conv_kelly_fft(test_mat = spatial_response_1, 
                                                     SF = SF_now)
      temporal_response_2 <- temporal_conv_kelly_fft(test_mat = spatial_response_2, 
                                                     SF = SF_now)
    }
    # checks:
    are_equal(temporal_response_1$size(), temporal_response_2$size())
    are_equal(temporal_response_1$size(), spatial_response_1$size())
    are_equal(temporal_response_2$size(), spatial_response_2$size())
    if (debug_mode) {
      print("Performed convolution through time, resulting in matrix of dim:")
      print(temporal_response_1$size())
      # progress of response through time
      temporal_resp_plots <- vector(mode = "list", length = 16)
      for (t_index_now_i in 1:16) {
        t_index_now <- seq(1, temporal_response_1$size(3), length.out = 16)[t_index_now_i]
        temporal_resp_plots[[t_index_now_i]] <- 
          plot_heatmap(as_array(temporal_response_1[ , , t_index_now]$cpu()), 
                       do_print = FALSE, do_rasterize = TRUE,
                       use_limits = c(as.numeric(torch_min(temporal_response_1)$cpu()), 
                                      as.numeric(torch_max(temporal_response_1)$cpu()))) + 
          theme(legend.position = "none") +
          ggtitle(paste("t =", round(t_index_now, 2) ))
      }
      temporal_resp_plots <- plot_grid(plotlist = temporal_resp_plots, nrow = 4, align = "hv")
      print(temporal_resp_plots)
      # # test for a single pixel:
      # plot(mat_size_t, as_array(spatial_response_2[80, 100, ]$cpu()), type = "l")
      # lines(mat_size_t, as_array(temporal_response_2[80, 100, ]$cpu()), col = "red")
    } else {
      # save memory:
      if (!is.null(irf_df)) {
        rm(resp_fun_gpu)
      }
      rm(spatial_response_1, spatial_response_2) 
    }
    ## 5. Squaring of zero-phase and half-phase responses
    square_override <- FALSE # squaring override? This should be FALSE
    if (!square_override) {
      squared_response <- torch_sqrt(torch_square(temporal_response_1) + torch_square(temporal_response_2))
    } else {
      squared_response <- torch_sqrt(torch_square(temporal_response_2))
    }
    if (debug_mode) {
      squared_resp_plots <- vector(mode = "list", length = 16)
      for (t_index_now_i in 1:16) {
        t_index_now <- seq(1, squared_response$size(3), length.out = 16)[t_index_now_i]
        squared_resp_plots[[t_index_now_i]] <- 
          plot_heatmap(as_array(squared_response[ , , t_index_now]$cpu()), 
                       do_print = FALSE, 
                       use_limits = c(0, as.numeric(torch_max(squared_response)$cpu()))) + 
          theme(legend.position = "none") +
          ggtitle(paste("t =", round(t_index_now, 2) ))
      }
      squared_resp_plots <- plot_grid(plotlist = squared_resp_plots, nrow = 4, align = "hv")
      print(squared_resp_plots)
      # max value?
      max(as.vector(as_array(squared_response[ , , t_index_now]$cpu())))
    } else {
      rm(temporal_response_1, temporal_response_2)
    }
    ## 6. prepare Delayed Normalization - LOWPASS Filter
    ## The squared response (all positive) will be convoluted with a declining-exponential
    ## low-pass filter. Then, this signal will serve as denominator
    # squared_response is now [y, x, t]
    # get the low-pass filter response function and reshape accordingly
    lp_fun_gpu <- torch_tensor(data = IRF_lp, device = use_device, dtype = preferred_dtype)
    lp_fun_gpu <- lp_fun_gpu$expand(c(squared_response$size(2), 1, lp_fun_gpu$size())) # expand to match matrix dimensions
    lp_fun_gpu <- torch_cat(list(torch_zeros(c(lp_fun_gpu$size(1), lp_fun_gpu$size(2), lp_fun_gpu$size(3)-1), 
                                             device = use_device, dtype = preferred_dtype), # create the zero pad to center the TRF
                                 lp_fun_gpu), dim = 3)
    if (debug_mode) {
      print("Created low-pass temporal response tensor of dimensions for delayed normalization:")
      print(lp_fun_gpu$size())
    }
    # convolve squared responses with low-pass filter
    lowpass_response <- torch_conv1d(input = squared_response, 
                                     weight = lp_fun_gpu$flip(3), 
                                     padding = lp_fun_gpu$size(3)%/%2, 
                                     groups = squared_response$size(2) # is the x dimension
    )
    are_equal(dim(lowpass_response), dim(squared_response))
    ## 6.1 add SMOOTHING IN SPACE to achieve a spatially distributed normalization pool
    if (use_normalization_pool) {
      ## first, create the smoothing kernel
      normalization_kernel_sd <- 7.5 # in pixels, roughly 0.5 dva if scr.ppd=15 
      normalization_kernel_size <- seq(round(-3*normalization_kernel_sd), 
                                       round(3*normalization_kernel_sd), 
                                       spatial_resolution) # this will always yield odd dimensions
      meshlist <- meshgrid(x = normalization_kernel_size, y = normalization_kernel_size)
      normalization_kernel <- exp(-(meshlist$X^2+meshlist$Y^2)/(2*normalization_kernel_sd^2)) 
      normalization_kernel <- normalization_kernel / sum(normalization_kernel)
      # second, convolve the low-pass filtered response with that kernel (which needs some moving dimensions again)
      lowpass_response <- lowpass_response$permute(c(3, 1, 2))$unsqueeze(1) # new dim: 1, t, y, x
      number_of_spatial_operations_normalization <- prod(dim(normalization_kernel)) * lowpass_response$size(3) * lowpass_response$size(4)
      lowpass_response_s <- NULL
      ## try with fft2 convolution, if it is favorable
      if ( (cuda_is_available() && number_of_spatial_operations_normalization > 10^9 ) ||
           ((!cuda_is_available() || no_CUDA) && number_of_spatial_operations_normalization > 10^7 ) 
      ) {
        tryCatch({
          lowpass_response_s <- fft_conv_with_torch_3d(x = lowpass_response, K = normalization_kernel, use_torch = TRUE,
                                                       no_CUDA = no_CUDA, use_half_precision = use_half_precision)
        }, error = function(e) { print(substr(e$message, 1, 250)) } )
      } 
      # in case this did not work, run the normal convolution
      if (is.null(lowpass_response_s)) {
        lowpass_response_s <- standard_conv_with_torch_3d(x = lowpass_response, K = normalization_kernel, 
                                                          no_CUDA = no_CUDA, use_half_precision = use_half_precision)
      }
      lowpass_response <- lowpass_response_s # overwrite lowpass response with spatially filtered one
      rm(lowpass_response_s)
      lowpass_response <- lowpass_response$permute(c(1, 3, 4, 2))$squeeze() # new dim: y, x, t
      are_equal(dim(lowpass_response), dim(squared_response)) # make sure we're back to the initial array dimension
      # clean up 
      if (!debug_mode) {
        rm(normalization_kernel_sd, normalization_kernel_size, meshlist, normalization_kernel, 
           number_of_spatial_operations_normalization)
      }
    }
    # have a look at the lowpass-filtered response:
    if (debug_mode) {
      array_lowpass_response <- as_array(lowpass_response$cpu())
      array_squared_response <- as_array(squared_response$cpu())
      where_max <- which(array_lowpass_response==max(array_lowpass_response), arr.ind=TRUE)[1, ]
      par(mfrow = c(1, 1))
      plot(mat_size_t, array_squared_response[where_max[1], where_max[2], ], type = "l", 
           ylim = c(0, max(array_squared_response)))
      lines(mat_size_t, array_lowpass_response[where_max[1], where_max[2], ], col = "blue")
      par(mfrow = c(1, 1))
      rm(array_lowpass_response, array_squared_response, where_max)
    } else {
      rm(lp_fun_gpu)
    }
    ## 7. save squared responses and lowpass-filtered responses 
    # we'll keep those in torch format, but store them in working memory where we have more space.
    # we'll store as half precision to save space...
    # and we'll also save a lot of computational time later, if we can just move it to GPU
    store_as <- torch_float16() # or torch_float16() if you want to save memory
    squared_3d_over_filters[[gabor_list_i]] <- squared_response$to(device = 'cpu', dtype = store_as)
    lowpass_3d_over_filters[[gabor_list_i]] <- lowpass_response$to(device = 'cpu', dtype = store_as)
    SFs_of_filters[gabor_list_i] <- as.numeric(SF_now)
    Oris_of_filters[gabor_list_i] <- as.numeric(Ori_now)
    toc(quiet = !(debug_mode || show_final_maps)) # keep track of time
    if (!debug_mode) {
      rm(squared_response, lowpass_response)
    } 
    
  } # END OF LOOP OVER GABOR FILTERS
  
  ## 8. NORMALIZATION
  if (debug_mode || show_final_maps) {
    if (!normalize_override) {
      if (use_normalization_pool) {
        print("Performing normalization using normalization pool...")
      } else {
        print("Performing normalization...")
      }
    }
  }
  # determine overall largest response
  max_squared_resp <- unlist(lapply(X = squared_3d_over_filters, 
                                    FUN = function(x) { as.numeric(max(x)) }))
  if (debug_mode) {
    max_squared_df <- data.table(SFs = SFs_of_filters, Oris = Oris_of_filters, 
                                 max_squared_resp = max_squared_resp)
    p_max_squared_df <- ggplot(data = max_squared_df, aes(x = Oris, y = SFs, 
                                                          fill = max_squared_resp)) + 
      geom_tile() + 
      coord_cartesian(expand = FALSE) + scale_y_log10() + 
      scale_fill_viridis_c()
    p_max_squared_df
  }
  # now determine the sigma for normalization: We'll scale it according to the maximum value in the set
  prm.sigma = 0.07 # again, see Figure 11 from "Temporal Dynamics of Neural Responses in Human Visual Cortex"
  prm.n = 1.4
  which_max_squared_resp <- which.max(max_squared_resp)
  if (is.null(norm_sigma)) { # no norm_sigma provided, so go estimate
    norm_sigma <- torch_scalar_tensor(value = max_squared_resp[which_max_squared_resp] * prm.sigma, 
                                      device = use_device, dtype = preferred_dtype)
  } else { # norm_sigma provided
    norm_sigma <- torch_scalar_tensor(value = norm_sigma, 
                                      device = use_device, dtype = preferred_dtype)
  }
  # preallocate the final output matrix [filter, y, x]
  final_2d_over_filters <- array(data = NaN, 
                                 dim = c(length(gabor_list), dim(mat_over_t)[1], dim(mat_over_t)[2]))
  # normalized final output huge matrix [filter, y, x, t]
  if (output_full_sequences) {
    required_dims <- c(length(gabor_list), dim(mat_over_t)[1], dim(mat_over_t)[2], dim(mat_over_t)[3])
    required_dimnames <- list(paste(SFs_of_filters, Oris_of_filters, sep = "_"),
                              mat_size_y, mat_size_x, mat_size_t)
    if (output_full_sequences_use_RleArray) { # use Run Length Encoding, sparse in memory, extremely recommended
      if ((debug_mode || show_final_maps)) { print("Creating RleArray with dimensions:") }
      dat <- lapply(X = 1:required_dims[1], 
                    FUN = function(x) { Rle(values = NaN, lengths = prod(required_dims[2:4])) })
      final_3d_over_filters <- RleArray(data = dat, 
                                        dim = required_dims, dimnames = required_dimnames)
      if ((debug_mode || show_final_maps)) { print(dim(final_3d_over_filters)) }
      # in case we have a delayed array, we must compress from time to time
      compress_how_often <- 0
      if (compress_how_often >= 2) {
        compress_when <- as.vector(round(seq(1, length(gabor_list), length.out = compress_how_often))[-1])
      } else if (compress_how_often == 1) {
        compress_when <- as.vector(length(gabor_list))
      } else if (compress_how_often == 0) {
        compress_when <- NULL
      }
      rm(dat)
    } else { # normal array: massive amount of memory needed
      final_3d_over_filters <- array(data = NaN,
                                     dim = required_dims, dimnames = required_dimnames )
    }
    rm(required_dims, required_dimnames)
  } else {
    final_3d_over_filters <- NULL
  }
  # for each single filter we will perform normalization based on the above values
  for (gabor_list_i in 1:length(gabor_list)) {
    tic(msg = paste0("Normalization of responses (filter ", gabor_list_i, " of ", length(gabor_list), ")"))
    ## Normalization here using torch and the GPU:
    if (!normalize_override) {
      # choose the normalization (using the lowpass-filtered data) pool here:
      if (use_normalization_pool) { # use normalization pool described in Schuett & Wichmann
        # parameters
        use_device_normalization_pool <- use_device #'cpu' # as transferring to gpu iteratively is costly
        n_adjacents <- 2 # how many adjacent values around the current
        sigma_SF <- 1 # octaves
        sigma_Ori <- 20 # degrees
        # what is the current SF and Orientation of the filter?
        norm_SF_now <- SFs_of_filters[gabor_list_i]
        unique_SFs <- sort(unique(SFs_of_filters))
        norm_Ori_now <- Oris_of_filters[gabor_list_i]
        unique_Oris <- sort(unique(Oris_of_filters))
        # get adjacent filter fields in frequency.
        adjacent_SFs_index <- (which(unique_SFs==norm_SF_now)-n_adjacents):
          (which(unique_SFs==norm_SF_now)+n_adjacents)
        adjacent_SFs_index <- adjacent_SFs_index[adjacent_SFs_index>=1 & 
                                                   adjacent_SFs_index<=length(unique_SFs)]
        adjacent_SFs <- unique_SFs[adjacent_SFs_index]
        freq_dist <- log2(adjacent_SFs)-log2(norm_SF_now)
        # get adjacent filter fields in orientation. 
        # this procedure assumes that the definition of orientations is truly circular, 
        # such as: -45, 0, 45, 90
        adjacent_Oris_index <- 
          (which(unique_Oris==norm_Ori_now)-n_adjacents):(which(unique_Oris==norm_Ori_now)+n_adjacents)
        adjacent_Oris_index[adjacent_Oris_index<1] <- 
          (length(unique_Oris)-sum(adjacent_Oris_index<1)+1):length(unique_Oris)
        adjacent_Oris_index[adjacent_Oris_index>length(unique_Oris)] <- 
          1:(sum(adjacent_Oris_index>length(unique_Oris)))
        adjacent_Oris <- unique_Oris[adjacent_Oris_index]
        ori_dist = norm_Ori_now - adjacent_Oris # 
        ori_dist[ori_dist<=(-pi/2)] <- ori_dist[ori_dist<=(-pi/2)] + pi
        ori_dist[ori_dist>=(+pi/2)] <- ori_dist[ori_dist>=(+pi/2)] - pi
        # determine weights based on the Gaussian distributions specified above
        SF_weights <- dnorm(x = freq_dist, mean = 0, sd = sigma_SF)
        SF_weights <- SF_weights / sum(SF_weights)
        Ori_weights <- dnorm(x = ori_dist, mean = 0, sd = deg2rad(sigma_Ori))
        Ori_weights <- Ori_weights / sum(Ori_weights)
        # get list indices
        norm_index <- which(is.element(SFs_of_filters, adjacent_SFs) & is.element(Oris_of_filters, adjacent_Oris))
        # to check: rbind(SFs_of_filters[norm_index], Oris_of_filters[norm_index])
        # ... and now walk through them
        normalize_by <- NULL
        for (index_now in norm_index) {
          SF_weight_now <- SF_weights[adjacent_SFs==SFs_of_filters[index_now]]
          Ori_weight_now <- Ori_weights[adjacent_Oris==Oris_of_filters[index_now]]
          assert_that(length(SF_weight_now)==1 & length(Ori_weight_now)==1)
          if (index_now==norm_index[1]) { # first set of data to enter normalization
            normalize_by <- lowpass_3d_over_filters[[index_now]]$to(device = use_device_normalization_pool, 
                                                                    dtype = preferred_dtype)  
            normalize_by <- torch_multiply(normalize_by, SF_weight_now*Ori_weight_now)
          } else { # add next set of data 
            normalize_by <- torch_add(normalize_by, 
                                      torch_multiply(lowpass_3d_over_filters[[index_now]]$to(device = use_device_normalization_pool, 
                                                                                             dtype = preferred_dtype), 
                                                     SF_weight_now*Ori_weight_now))
          }
        }
        normalize_by <- normalize_by$to(device = use_device,
                                        dtype = preferred_dtype)
      } else { # NO normalization pool, simply the original delayed normalization model
        normalize_by <- lowpass_3d_over_filters[[gabor_list_i]]$to(device = use_device, 
                                                                   dtype = preferred_dtype)
      }
      # DIVISIVE NORMALIZATION:
      normalized_response <- 
        torch_divide(torch_pow(squared_3d_over_filters[[gabor_list_i]]$to(device = use_device, 
                                                                          dtype = preferred_dtype), 
                               prm.n), 
                     torch_pow(norm_sigma, prm.n) + torch_pow(normalize_by, prm.n) )
    } else { # NO NORMALIZATION, as normalize_override==TRUE
      normalized_response <- squared_3d_over_filters[[gabor_list_i]]$to(device = use_device, 
                                                                        dtype = preferred_dtype)
    }
    are_equal(dim(normalized_response), dim(lowpass_3d_over_filters[[gabor_list_i]]))
    are_equal(dim(squared_3d_over_filters[[gabor_list_i]]), dim(lowpass_3d_over_filters[[gabor_list_i]]))
    # save this response for plotting later
    if (#debug_mode && 
        gabor_list_i==which_max_squared_resp) {
      max_normalized_filter <- as_array(normalized_response$cpu())
    }
    ## Sum-over-time operation, reduce 3D to 2D [y, x]
    normalized_response_sum <- normalized_response$sum(dim = 3)
    ## Pack into final array
    final_2d_over_filters[gabor_list_i, , ] <- as_array(normalized_response_sum$cpu())
    if (output_full_sequences) {
      # this assignment is either delayed or not:
      final_3d_over_filters[gabor_list_i, , , ] <- as_array(normalized_response$cpu())
      # if it is delayed, then we must compress from time to time:
      if (output_full_sequences_use_RleArray && !is.null(compress_when) && (gabor_list_i %in% compress_when)) {
        if ((debug_mode || show_final_maps)) { print("Compressing final_3d_over_filters ...") }
        final_3d_over_filters <- as(final_3d_over_filters, "RleArray")
        if ((debug_mode || show_final_maps)) { print("Done.") }
      }
    }
    # 
    if (debug_mode) {
      print(paste("Computed sum. Resulting dimensions:", paste0(dim(normalized_response_sum), collapse = ",")))
      p_sum <- plot_heatmap(mat = final_2d_over_filters[gabor_list_i, , ], do_print = FALSE)
      print(p_sum)
    } else {
      rm(normalized_response_sum, normalized_response)
      if (!normalize_override) { 
        rm(normalize_by) 
        if (use_normalization_pool) {
          rm(adjacent_SFs_index, adjacent_Oris_index, ori_dist, freq_dist, 
             SF_weights, Ori_weights, norm_index)
        }
      }
    }
    # take time
    toc(quiet = !(debug_mode || show_final_maps))
  } # end of loop over filters, for which normalization is performed
  
  # diagnostic plots on normalization?
  if (debug_mode) {
    # get squared, lowpass filters and find out the y,x position with max response
    max_squared_filter <- as_array(squared_3d_over_filters[[which_max_squared_resp]]$cpu())
    max_lowpass_filter <- as_array(lowpass_3d_over_filters[[which_max_squared_resp]]$cpu())
    where_max <- which(max_squared_filter==max(max_squared_filter), arr.ind=TRUE)
    # normalize the normalized data, so that it fits the same plot
    norm_max_normalized_filter <- max_normalized_filter[where_max[1], where_max[2], ]
    norm_max_normalized_filter <- norm_max_normalized_filter / max(norm_max_normalized_filter) * max(max_squared_filter)
    # plot now:
    plot(mat_size_t, max_squared_filter[where_max[1], where_max[2], ], type = "l")
    lines(mat_size_t, max_lowpass_filter[where_max[1], where_max[2], ], col = "blue")
    lines(mat_size_t, norm_max_normalized_filter, col = "red")
    rm(max_squared_filter, max_lowpass_filter, max_normalized_filter, 
       norm_max_normalized_filter, where_max)
  } else {
    rm(squared_3d_over_filters, lowpass_3d_over_filters)
  }
  # clean up at the end
  if (!debug_mode) {
    rm(mat_over_t, mat_over_t_gpu)
  }
  # a final matrix of 2D responses?
  if (show_final_maps && final_maps_filename!="") {
    print("Plotting...")
    all_final_maps <- vector(mode = "list", length = length(gabor_list))
    for (gabor_list_i in 1:length(gabor_list)) {
      all_final_maps[[gabor_list_i]] <- plot_heatmap(mat = final_2d_over_filters[gabor_list_i, , ], 
                                                     do_print = FALSE, do_rasterize = FALSE,
                                                     use_grayscale = FALSE, use_viridis = TRUE,
                                                     use_limits = c(0, max(final_2d_over_filters))) + 
        ggtitle(paste0("SF=", round(SFs_of_filters[gabor_list_i], 2), 
                       ", ",
                       "Ori=", round(Oris_of_filters[gabor_list_i], 2) )) + 
        theme_void() + theme(legend.position = "none")
    }
    require(gridExtra)
    ggsave(final_maps_filename, arrangeGrob(grobs = all_final_maps, 
                                            nrow = length(unique(SFs_of_filters)), 
                                            ncol = length(unique(Oris_of_filters))), 
           device = "pdf", 
           width = 15, height = 12)
  } else {
    all_final_maps <- NULL
  }
  if (debug_mode || show_final_maps) { 
    print("Done.")
  }
  # return
  return(list(final_2d_over_filters, SFs_of_filters, Oris_of_filters, final_3d_over_filters, all_final_maps))
}
