multinom_sampler <- function(cat_prob, 
                             random_seed = NaN, 
                             n_samples = 100000,
                             categories = c("static", "inward", "outward", "downward", "upward"), 
                             do_print = FALSE) {
  # if the supplied category probabilities should exceed 1, normalize to 1
  if (sum(cat_prob)>1) { 
    cat_x <- categories
    cat_prob_x <- cat_prob / sum(cat_prob)
  } else { # there's room for guessing
    cat_x <- c(categories, "nowletusguess")
    cat_prob_x <- c(cat_prob, 1-sum(cat_prob))
  }
  # get multinomial sample
  if (!is.na(random_seed)) { set.seed(1000+random_seed) }
  sam <- sample(cat_x, n_samples, 
                replace = TRUE, 
                prob = cat_prob_x)
  # for those where let us guess occurs, sample again
  if (!is.na(random_seed)) { set.seed(1000000+random_seed) }
  sam[sam=="nowletusguess"] <- sample(categories, sum(sam=="nowletusguess"), 
                                      replace = TRUE, 
                                      prob = rep(1, length(categories)) / length(categories))
  if (do_print) {
    print(rbind(categories, cat_prob))
    print(prop.table(table(sam, useNA = "ifany")))
  }
  return(sam)
}
