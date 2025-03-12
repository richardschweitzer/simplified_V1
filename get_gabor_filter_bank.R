get_gabor_filter_bank <- function(all_SF, all_Ori, scr_ppd, ppd_scaler = 1, 
                                  use_log_gabor_heiko = FALSE, # use Heiko Schuett's version?
                                  use_bw_heiko = c(0.5945, 0.2965), # SDs for SF and Orientation in log frequency space, according to Heiko
                                  use_log_gabor_SF_band_adjustment = FALSE, # if log Gabor, then SF bandwidth adjustment at low SFs?
                                  use_log_gabor_orientation_band_adjustment = FALSE, # if log Gabor, then orientation bandwidth adjustment?
                                  do_on_GPU=FALSE, # transfer Gabor filters to GPU memory?
                                  rf_width_override = NaN,
                                  make_gaussian_aperture = FALSE, 
                                  do_plot=FALSE) {
  if (do_on_GPU) {
    require(torch) # we can use torch to save these gabors on the GPU memory
    if (cuda_is_available()) {
      devi <- 'cuda'
    } else {
      devi <- 'cpu'
    }
  }
  # CREATE A BANK OF GABOR FILTERS
  resulting_gabors <- vector(mode = "list", length = length(all_SF)*length(all_Ori))
  result_i <- 0
  for (now_SF_i in 1:length(all_SF) ) {
    now_SF <- all_SF[now_SF_i]
    # have we specified the RF size?
    if (any(is.na(rf_width_override)) || length(rf_width_override)==1) {
      rf_width_override_now <- rf_width_override
    } else if (length(rf_width_override)==length(all_SF)) {
      rf_width_override_now <- rf_width_override[now_SF_i]
    }
    for (now_Ori in all_Ori) {
      result_i <- result_i + 1
      if (use_log_gabor_heiko) {
        # perform adjustment of SF bandwidth? 
        if (use_log_gabor_SF_band_adjustment) {
          assert_that(exists("anderson_burr_SF_bandwidth"))
          assert_that(exists("SD_to_bandwidth_at_half_height"))
          # 1. convert Heiko's Gaussian SD to full bandwidth
          minimum_bandwidth_now <- SD_to_bandwidth_at_half_height(use_bw_heiko[1])
          # 2. determine bandwidth of current SF
          bandwidth_now = anderson_burr_SF_bandwidth(SF = now_SF, 
                                                     bw_minimum = minimum_bandwidth_now)
          # 3. transform full bandwidth back to Gaussian SD
          use_SF_bw = SD_to_bandwidth_at_half_height(bandwidth_now, run_inverse = TRUE)
        } else {
          use_SF_bw = use_bw_heiko[1]
        }
        # perform adjustment of orientation bandwidth? 
        if (use_log_gabor_orientation_band_adjustment) {
          assert_that(exists("anderson_burr_RF_width"))
          assert_that(exists("orientation_bandwidth_from_aperture_SD"))
          assert_that(exists("SD_to_bandwidth_at_half_height"))
          # 1. get theoretical RF width
          RF_width_now <- anderson_burr_RF_width(now_SF) 
          # 2. get theoretical bandwidth based on Gabor filter definition: SF and Gaussian envelope SD
          RF_bandwidth_now <- orientation_bandwidth_from_aperture_SD(SF = now_SF, 
                                                                     Gaussian_SD = RF_width_now/2)
          # 3. convert full bandwidth to Gaussian SD
          use_Ori_bw = SD_to_bandwidth_at_half_height(RF_bandwidth_now, run_inverse = TRUE)
          #print(paste(round(now_SF, 2), round(rad2deg(RF_bandwidth_now), 2), round(use_Ori_bw, 3)))
        } else {
          use_Ori_bw = use_bw_heiko[2]
        }
        # create log Gabor with specified bandwidths
        this_gabor <- get_gabor_field_heiko(rf_freq_dva = now_SF, rf_ori = now_Ori, 
                                            bw = c(use_SF_bw, use_Ori_bw),
                                            im_size_dva = rf_width_override_now, 
                                            scr.ppd = scr_ppd, ppd_scaler = ppd_scaler)
      } else {
        this_gabor <- get_gabor_field(rf_freq_dva = now_SF, rf_ori = now_Ori, 
                                      create_aperture = TRUE, 
                                      gaussian_aperture = make_gaussian_aperture,
                                      rf_width_dva = rf_width_override_now, 
                                      scr.ppd = scr_ppd, ppd_scaler = ppd_scaler)
      }
      # do we want to perform this operation on the GPU? If so, transfer to GPU using torch
      if (do_on_GPU) {
        this_gabor[[1]] <- torch_tensor(this_gabor[[1]], device = devi)
        this_gabor[[2]] <- torch_tensor(this_gabor[[2]], device = devi)
      }
      resulting_gabors[[result_i]] <- this_gabor
      if (do_plot & !do_on_GPU) {
        p_gabor <- plot_heatmap(resulting_gabors[[result_i]][[1]])
        print(p_gabor)
      }
    }
  }
  return(resulting_gabors)
}