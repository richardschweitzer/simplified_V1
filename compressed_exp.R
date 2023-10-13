# classic compressed exp model without delay
compressed_exp <- function(t, delay=0, amp, dur, tail) { 
  res = amp * (1 - exp(-1 * ((t-delay)/dur)^tail))
  res[t<delay] <- 0
  return(res)
}