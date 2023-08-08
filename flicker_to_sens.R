flicker_to_sens <- function(tfreqs,  # can be numeric vector
                            a0, a1, a2, a3, a4 = NaN, # scalar values
                            # if a4 is set to NaN, then Burr & Morrone function will be used, 
                            # otherwise, Bergen & Wilson will be used, where a0=A, a1=B, a2=tau, a3=n, a4=k
                            t_res, pres_dur, ramp_sd=NaN, # those are constants
                            beta=4, G=1 # slope of prob summation, negative-lobe asymmetry
) {
  # create temporal reference
  stim_time <- seq(0, pres_dur, by = t_res) 
  # iterate over tfreqs
  sums_tfreqs <- vector(mode = "numeric", length = length(tfreqs))
  for (tfreq_i in (1:length(tfreqs))) {
    # create flicker stimulus and gaussian temporal ramp
    x <- create_flicker_stim(dur = pres_dur, t_res = t_res, 
                             tfreq = tfreqs[tfreq_i]) 
    if (!is.na(ramp_sd)) {
      gaussian <- dnorm(x = stim_time, mean = pres_dur/2, sd = ramp_sd) # gaussian ramp
      gaussian <- gaussian / max(gaussian)
    } else {
      gaussian <- rep(1, length.out = length(stim_time))
    }
    stim_tfreq <- x * gaussian # combined
    # create IRF
    if (is.na(a4)) {
      trf <- exp_damp_sin(x_step = t_res, x_dur = 500, latency = 0, 
                          a0 = a0, a1 = a1, a2 = a2, a3 = a3)
    } else {
      trf <- bergen_wilson(x_step = t_res, x_dur = 500, latency = 0, 
                           A = a0, B = a1, tau = a2, n = a3, k = a4)
    }
    trf_time <- trf[[2]]
    trf <- trf[[1]]
    # convolution
    resp_tfreq <- zapsmall(convolve(stim_tfreq, rev(trf), type = "open"))[1:length(stim_tfreq)]
    # add the systematic asymmetry? (as in Kelly & Savoie, 1978 or Bergen & Wilson, 1984)
    if (G != 1) {
      resp_tfreq[resp_tfreq<0] <- resp_tfreq[resp_tfreq<0] * G
    }
    # probability summation
    sums_tfreqs[tfreq_i] <- pracma::trapz(x = stim_time, y = abs(resp_tfreq)^beta)^(1/beta)
  }
  return(sums_tfreqs)
}
