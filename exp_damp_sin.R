# this is the function for luminance with parameters of DB
exp_damp_sin_fun <- function(t_ms, a0 = 1.8e5, a1 = 13, a2 = 5, a3 = 27.3) {
  return(
    a0 * as.numeric(t_ms/1000>=0) * t_ms/1000 * sin(2*pi*(a1*t_ms/1000*(t_ms/1000+1)^(-a2))) * exp(-a3*t_ms/1000) 
  ) 
}

exp_damp_sin <- function(x_dur = 300, x_step = 1, # in milliseconds
                         a0 = 1.8e5, a1 = 13, a2 = 5, a3 = 27.3, 
                         latency = 0, norm_amp = NaN) {
  # exponentially damped freq-modulated sinusoid (Burr & Morrone, 1993)
  x_here <- seq(from = 0, to = x_dur, by = x_step)
  res <- exp_damp_sin_fun(t_ms = x_here, a0 = a0, a1 = a1, a2 = a2, a3 = a3)
  if (latency>0) { # apply additional latency?
    latency_here <- seq(from = 0, to = latency, by = x_step)
    res <- c(rep(0, length(latency_here)), res)
    x_here <- c(latency_here, x_here+max(latency_here)+x_step)
  }
  if (!is.na(norm_amp)) { # normalize amplitude to one?
    res <- res / max(res)
    res <- res * norm_amp
  }
  return(list(res, x_here))
}
