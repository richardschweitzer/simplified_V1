


####################### prerequisites ############################
options(torch.threshold_call_gc = 6000)
library(torch)
library(assertthat)
library(data.table)
library(ggplot2)
library(cowplot)
library(scales)
library(viridis)
library(cowplot)
library(ggforce)
library(scales)
library(pracma)
library(mgcv)
library(tictoc)
library(minpack.lm)
library(WRS2)
library(lme4)
# deal with massive objects in memory, required: install.packages("BiocManager")
# then: BiocManager::install("DelayedArray")
# library(DelayedArray)
# here: BiocManager::install("HDF5Array")
# library(HDF5Array)
# or here: BiocManager::install("SparseArray")
# library(SparseArray)

SDECTheme <- function(base_size=15, base_family="Helvetica") { 
  theme_classic(base_size=base_size, base_family=base_family) %+replace% 
    theme(
      # size of text
      axis.text = element_text(size = base_size), # 0.9*base_size
      legend.text = element_text(size = base_size),
      # remove grid horizontal and vertical lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background  = element_blank(),
      # panel boxes
      panel.border = element_blank(),
      axis.line = element_line(colour = "grey20"), 
      # facet strips
      strip.text = element_text(size = base_size, face = "bold"),
      strip.background = element_rect(fill="transparent", colour = "transparent"),
      strip.placement = "outside",
      # legend
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA)
    )
}


# source the two spatial convolution functions here:
source("standard_conv_with_torch_3d.R")
source("fft_conv_with_torch_3d.R")
source("plot_heatmap.R")
source("temporal_conv_kelly_fft.R")
source("kelly_vel.R")
source("resample_yxt.R")

# this is the V1 module
source("v1.R")






##################### get an image #########################
require(jpeg)
library(colorspace)
source("plot_heatmap.R")

# read the image presented in the trial
stim_image <- readJPEG("example_image.jpg")
# convert to grayscale
stim_image_hex <- rgb(stim_image[,,1], stim_image[,,2], stim_image[,,3]) # hexadecimal color space
stim_image_hex_gray <- desaturate(stim_image_hex) # grayscale hexadecimal color space
stim_image_gray <- col2rgb(stim_image_hex_gray)[1, ]/255 # converted back to 0..1 colorspace
dim(stim_image_gray) <- dim(stim_image_hex) <- dim(stim_image_hex_gray) <- dim(stim_image)[1:2]
# trim the image to fullHD resolution
stim_image_pres_offset <- ceiling((dim(stim_image_gray) - c(1080, 1920)) / 2) # we will trim to this
# delete the old versions of the image
rm(stim_image, stim_image_hex, stim_image_hex_gray)
# pick the presented portion of the image
stim_image_gray <- stim_image_gray[stim_image_pres_offset[1]:(stim_image_pres_offset[1]+1080-1), 
                                   stim_image_pres_offset[2]:(stim_image_pres_offset[2]+1920-1)]
assert_that(all(dim(stim_image_gray)==c(1080, 1920)))
# show image
plot_heatmap(stim_image_gray, use_grayscale = TRUE)


#################### generate saccades ######################
# classic compressed exp model with delay
compressed_exp <- function(t, delay=0, amp, dur, tail) { 
  res = amp * (1 - exp(-1 * ((t-delay)/dur)^tail))
  res[t<delay] <- 0
  return(res)
}
# lebedev's main sequence function
lebedev_sqrt <- function(sac_amp, f1) {
  sac_dur <- f1 * sqrt(sac_amp)
  return(sac_dur)
}


# define amplitudes and durations, according to main sequence
amps <- seq(2, 20, by = 2)
durs <- lebedev_sqrt(amps, f1 = 15)

# time course
sac_pad = 150
ramp_dur = 50
time <- seq(-(ramp_dur+sac_pad), 50+sac_pad, by = 1)
detect_thres <- 0.01 # for velocity-based "saccade detection"
# create contrast ramp for onset. Linear?
time_contrast_ramp <- rep(1, length(time))
time_contrast_ramp[1:ramp_dur] <- linspace(0, 1, n = length(1:ramp_dur))
# create saccade trajectory
all_amp_traj <- vector(mode = "list", length = length(amps))
all_amp_onsets <- vector(mode = "numeric", length = length(amps))
all_amp_offsets <- vector(mode = "numeric", length = length(amps))
all_amp_vpeaks <- vector(mode = "numeric", length = length(amps))
for (amp_i in rev(1:length(amps))) {
  # create trajectory
  all_amp_traj[[amp_i]] <- compressed_exp(time, 
                                          amp = amps[amp_i], dur = durs[amp_i]/2, 
                                          tail = 2.5)
  # determine offset by threshold
  all_amp_onsets[amp_i] <- min(time[c(0, diff(all_amp_traj[[amp_i]])) > detect_thres], 
                               na.rm = TRUE)
  all_amp_offsets[amp_i] <- max(time[(c(0, diff(all_amp_traj[[amp_i]])) > detect_thres) &
                                       (time>all_amp_onsets[amp_i])], 
                               na.rm = TRUE)
  all_amp_vpeaks[amp_i] <- max(time[(c(0, diff(all_amp_traj[[amp_i]])) == max(diff(all_amp_traj[[amp_i]]))) &
                                       (time>all_amp_onsets[amp_i])], 
                                na.rm = TRUE)
  # plot
  if (amp_i==length(amps)) {
    plot(time, c(0, diff(all_amp_traj[[amp_i]])), type = "l", xlim = c(-20, 100), 
         ylab = "diff(traj)")
    abline(v = all_amp_offsets[amp_i])
  } else {
    lines(time, c(0, diff(all_amp_traj[[amp_i]])))
    abline(v = all_amp_offsets[amp_i])
  }
}

# summary
rbind(amps, durs, 
      all_amp_onsets, all_amp_offsets, all_amp_vpeaks)


#################### create the retinal input and run model ###########################
require(fields)
require(e1071) # for Hanning window

# parameters
scr_ppd <- 30
view_size = 10*scr_ppd # in pixels
(SFs = 1)#exp(linspace(log(0.25), log(4), 2))) # 7 originally, 0.25 - 4
(Oris = -pi/2)#seq(-pi/2, pi/2, length.out = 9)[-1]) # 9 originally

# create Gabor fields
source("get_log_gabor_heiko.R")
source("get_gabor_field_heiko.R")
source("get_gabor_field.R")
source("get_gabor_filter_bank.R")
source("anderson_burr_RF_width.R")
source("anderson_burr_SF_bandwidth.R")
source("SD_to_bandwidth_at_half_height.R") # conversion function Gaussian SD <> bandwidth at half height
source("orientation_bandwidth_from_aperture_SD.R")
# resulting from the underdetermination of the filters
gabor_list_heiko_2 <- get_gabor_filter_bank(
  use_log_gabor_heiko = TRUE,
  use_log_gabor_SF_band_adjustment = TRUE, # use SF adjustment (default: FALSE)
  use_log_gabor_orientation_band_adjustment = TRUE, # use Ori adjustment (default: FALSE)
  all_SF = SFs, 
  all_Ori = Oris, # positive means clockwise
  rf_width_override = anderson_burr_RF_width(SFs) * 8,
  scr_ppd = scr_ppd, 
  ppd_scaler = 1,
  make_gaussian_aperture = TRUE)
length(gabor_list_heiko_2)


# and TRFs
source("get_IRF_df.R")
irf_space <- get_IRF_df(SFs = SFs, irf_time_range = c(0, 200), 
                        temporal_resolution = 1)
ggplot(data = irf_space, aes(x = irf_time, y = irf)) + 
  geom_line(linewidth = 2) + 
  theme_minimal() + SDECTheme() + 
  facet_wrap(~round(SF, 2))


# make image object for interpolation
image_obj <- list(x = 1:dim(stim_image_gray)[1], # x and y must match dimension numbers?
                  y = 1:dim(stim_image_gray)[2], 
                  z = stim_image_gray) 
# create a Hanning window
hw <- hanning.window(view_size)
hw2 <- as_array(torch_multiply(t(t(hw)), t(hw)))
# interpolate view based on each retinal trajectory, then run model
model_output_over_time <- NULL
for (amp_i in rev(1:length(amps))) {
  print(paste("Running model for amplitude", amps[amp_i], "..."))
  # place trajectory around center of image
  eye_x = dim(stim_image_gray)[2]/4 +#- amps[amp_i]/2*scr_ppd + 
    all_amp_traj[[amp_i]] * scr_ppd
  eye_y = rep(dim(stim_image_gray)[1]/2, length(eye_x))
  # create sequence
  eye_x_seq <- seq(-view_size/2, view_size/2, length.out = view_size)
  eye_y_seq <- seq(-view_size/2, view_size/2, length.out = view_size)
  # now interpolate
  all_retinal_im <- array(data = NA, dim = c(view_size, view_size, length(eye_x)), 
                          dimnames = list(eye_y_seq, eye_x_seq, time)) # this array should have dimnames
  for (t_i in 1:length(eye_x)) {
    # get all points
    eye_x_points <- eye_x[t_i] + eye_x_seq
    eye_y_points <- eye_y[t_i] + eye_y_seq
    # expand
    query_points <- expand.grid(rows = eye_y_points, 
                                cols = eye_x_points)
    # interpolate
    retinal_im <- interp.surface( obj = image_obj, loc = query_points )
    retinal_im <- matrix(data = retinal_im, 
                         nrow = view_size, ncol = view_size, 
                         byrow = FALSE)
    # save:
    all_retinal_im[ , , t_i] <- retinal_im * hw2 * time_contrast_ramp[t_i] # and apply the Hanning window
  }
  # have a look
  if (FALSE) {
    for (plot_i in (c(290, 300, 310, 320, 330, 340, 350, 360, 370, 380)-150)) {
      plot_heatmap(all_retinal_im[ , , plot_i], use_grayscale = TRUE, do_print = TRUE)
    }
  }
  rm(retinal_im, query_points, eye_x_points, eye_y_points)
  
  
  # now we can run the model on this simulated retinal input
  test_RF_output <- v1(mat_over_t = all_retinal_im,
                       gabor_list = gabor_list_heiko_2, # a list of gabor filters created by 'get_gabor_filter_bank'
                       override_spatial_use_fft = FALSE, # to not use fft-based convolutions -> memory issue
                       irf_df = irf_space, # a data.frame of IRFs resolved in SF x time created by 'get_IRF_df'
                       normalize_override = FALSE,
                       use_normalization_pool = FALSE, 
                       output_full_sequences = TRUE, # TRUE: returns complete sequences. WARNING: this is memory-intense
                       output_full_sequences_spatial_resample_to = c(30, 30), 
                       norm_sigma = 0.07, 
                       signal_x_range = NULL, signal_y_range = NULL, signal_t_range = NULL, # the min and max of temporal scale
                       no_CUDA = FALSE, # set to TRUE if you want to use CPU instead of GPU (but torch)
                       use_half_precision = FALSE, # if a GPU is available, use half-precision?
                       spatial_resolution = 1, # spatial resolution of processing function
                       temporal_resolution = 1, # temporal resolution of processing function in milliseconds
                       debug_mode = FALSE, # shows output at every step
                       show_final_maps = TRUE,  # shows all resulting 2D maps at the end
                       final_maps_filename = "" # name of the pdf file if show_final_maps==TRUE
  )
  # get output and compress across channel and space
  out <- test_RF_output[[4]]
  out_over_channels <- apply(out, MARGIN = c(2, 3, 4), FUN = sum)
  out_over_time <- apply(out, MARGIN = c(4), FUN = sum)
  if (amp_i==length(amps)) {
    plot(time, out_over_time, type = "l")
  } else {
    lines(time, out_over_time)
  }
  
  # save
  model_output_over_time <- rbind(model_output_over_time, 
                                  data.table(time = time, 
                                             sac_on = all_amp_onsets[amp_i], 
                                             sac_off = all_amp_offsets[amp_i],
                                             sac_vpeak = all_amp_vpeaks[amp_i],
                                             eye_x = eye_x / scr_ppd, 
                                             eye_y = eye_y / scr_ppd, 
                                             amp = amps[amp_i],
                                             out = out_over_time))
  # clean up
  rm(test_RF_output, out, out_over_channels)
}

# final plot, saccade onset
ggplot(data = model_output_over_time, 
       aes(x = time, y = out, group = amp, color = amp)) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_line(size = 1.5, alpha = 0.8) + 
  SDECTheme() + 
  scale_color_viridis_c() + 
  labs(x = "Time re sac onset", y = "Model output sum", 
       color = "Sacc. amp")

# saccade vpeak
ggplot(data = model_output_over_time, 
       aes(x = time-sac_vpeak, y = out, group = amp, color = amp)) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_line(size = 1.5, alpha = 0.8) + 
  SDECTheme() + 
  scale_color_viridis_c() + 
  labs(x = "Time re sac peak velocity", y = "Model output sum", 
       color = "Sacc. amp")

# saccade offset
ggplot(data = model_output_over_time, 
       aes(x = time-sac_off, y = out, group = amp, color = amp)) + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_line(size = 1.5, alpha = 0.8) + 
  SDECTheme() + 
  scale_color_viridis_c() + 
  labs(x = "Time re sac offset", y = "Model output sum", 
       color = "Sacc. amp")





