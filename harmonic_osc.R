# harmonic oscillator step function (by Specht et al., 2017)
harmonic_osc <- function(t, delay=0, amp, dur, omega) {
  Xstep <- amp * (1 - exp(-1 * ((t-delay)/dur)) * (cos(omega*(t-delay)) + ((1/dur)/omega) * sin(omega*(t-delay))))
  Xstep[t<delay] <- 0
  return(Xstep)
}