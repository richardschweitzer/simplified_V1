get_IRF_df <- function(SFs, irf_time_range=c(0, 300), temporal_resolution, 
                       path_to=".", irf_filename="gam_irf.rda", 
                       make_sure_that_odd=FALSE) {
  require(mgcv)
  require(data.table)
  # load the TRF GAM
  load(file.path(path_to, irf_filename))
  # create time vector, make sure it has odd dimensions
  irf_time_now <- seq(from = irf_time_range[1], to = irf_time_range[2], by = temporal_resolution)
  if (make_sure_that_odd==TRUE && mod(length(irf_time_now), 2) == 0) {
    irf_time_now <- irf_time_now[-length(irf_time_now)]
  }
  # create space for the GAM prediction
  irf_space <- data.table(expand.grid(list(
    irf_time = irf_time_now, 
    SF = SFs
  )))
  # let the TRF GAM predict
  irf_space[ , irf := predict(gam_irf, newdata = irf_space)]
  # last checks
  assert_that(all(!is.na(irf_space$irf) & !is.infinite(irf_space$irf)))
  # output
  return(irf_space)
}