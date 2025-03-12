anderson_burr_SF_bandwidth <- function(SF, bw_minimum, SF_angle_point=2, do_plot=FALSE) {
  # bw_minimum: value of the bandwidth above the angle point 2 cpd
  # (Anderson & Burr have 1.5 octaves full bandwidth at half height)
  
  # to debug and to replicate Anderson's & Burr's results: 
  # SF = exp(linspace(log(0.03), log(30), 100))
  # bw_minimum = 1.5
  # to assert the correctness of these results, you should arrive at a bandwidth of 3 for 0.03 cpd
  
  # convert SF to octave scale relative to the angle point, i.e., 2 cpd
  SF_ap <- SF_angle_point
  oct <- log2(SF/SF_ap)
  # compute linear increase of "0.26 octaves per octave of spatial frequency"
  bandwidths <- rep(NaN, length(SF))
  bandwidths[oct>=0] <- bw_minimum
  bandwidths[oct<0] <- -0.26 * oct[oct<0] + bw_minimum
  # plot?
  if (do_plot) {
    require(ggplot2)
    df <- data.frame(SF = SF, oct = oct, bandwidths = bandwidths)
    p_df <- ggplot(df, aes(x = SF, y = bandwidths)) + 
      geom_line() + 
      scale_y_continuous(limits = c(0, 4)) + 
      scale_x_log10(limits = c(0.02, 50)) + 
      annotation_logticks(sides = "b") +
      theme_minimal()
    print(p_df)
  }
  # that's it
  return(bandwidths)
}