bergen_wilson_fun <- function(t_ms, A=21e3, B=0.92, tau=9.5, n=4, k=1.5) {
  # H(t) from Bergen & Wilson (1984), with standard params for sigma=1
  return(
    A * (t_ms/tau)^n * exp(-(t_ms/tau)) * (1/factorial(n) - B * ((t_ms/tau)^k) / factorial(n+k) )
  )
}

bergen_wilson <- function(x_dur = 300, x_step = 1, # in milliseconds
                          A = 21e3, B = 0.92, tau = 9.5, n = 4, k = 1.5,
                          latency = 0, norm_amp = NaN) {
  # exponentially damped freq-modulated sinusoid (Burr & Morrone, 1993)
  x_here <- seq(from = 0, to = x_dur, by = x_step)
  res <- bergen_wilson_fun(t_ms = x_here, A = A, B = B, tau = tau, n = n, k = k)
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
