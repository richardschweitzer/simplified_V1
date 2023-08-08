# function to create flicker stimulus
create_flicker_stim <- function(dur, t_res, tfreq) {
  T_ms <- 1000 / tfreq
  soa_tfreq <- T_ms / 2
  stim_time <- seq(0, dur, by = t_res)
  stim_polarity <- sin(2*pi*stim_time/1000*tfreq)
  return(stim_polarity)
}
