load_noise_patch <- function(path_to="stimuli_matrices", stim_name, 
                             subtract_to_zero=0.5, do_clip=0.5, upscale_fac=2) {
  noise_patch <- as.matrix(read.csv(file.path("stimuli_matrices", stim_name), header = FALSE))
  (noise_patch_dim <- dim(noise_patch))
  dimnames(noise_patch) <- NULL
  if (!is.na(subtract_to_zero)) {
    noise_patch <- noise_patch - subtract_to_zero
    if (!is.na(do_clip)) {
      noise_patch[noise_patch>do_clip] <- do_clip
      noise_patch[noise_patch<(-1*do_clip)] <- (-1)*do_clip
    }
    if (!is.na(upscale_fac)) {
      noise_patch <- noise_patch * upscale_fac
    }
  }
  return(noise_patch)
}
