plot_heatmap <- function(mat, 
                         use_grayscale=FALSE, 
                         use_viridis=FALSE,
                         use_limits=NULL, 
                         do_print=FALSE, 
                         do_rasterize=FALSE) { # richards function: plots an image matrix with ggplot2
  require(reshape2)
  require(ggplot2)
  if (use_viridis) {
    require(viridis)
  }
  if (do_rasterize) {
    require(ggrastr)
  }
  mat_long <- reshape2::melt(mat)
  colnames(mat_long) <- c("y", "x", "z")
  p <- ggplot(mat_long, aes(x = x, y = y, fill = z))
  if (do_rasterize) {
    p <- p + rasterise(geom_tile(), dpi = 96) 
  } else {
    p <- p + geom_tile() 
  }
  p <- p + theme_minimal() + 
    coord_fixed(expand = FALSE) + scale_y_reverse()
  if (use_grayscale) {
    p <- p + scale_fill_distiller(type = "seq",
                                  direction = -1,
                                  palette = "Greys", 
                                  limits = use_limits)
  } else {
    if (use_viridis) {
      p <- p + scale_fill_viridis_c(limits = use_limits)
    } else {
      p <- p + scale_fill_gradient2(limits = use_limits)
    }
  }
  if (do_print) {
    print(p)
  }
  return(p)
}