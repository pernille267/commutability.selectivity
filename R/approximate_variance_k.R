#' Approximate variance of k given number of clinical samples and replicates
#'
#' @param n Integer - The number of unique clinical samples
#' @param R Integer - The number of replicated measurements done on each clinical sample
#' @param M Double - The relative increase in prediction interval width between the case of identical selectivity profiles and any other case
#'
#' @return A value that approximates the true variance of k
#' @export
#'
#' @examples approximate_variance_k(25, 3)

approximate_variance_k <- function(n, R, M = 0){
  expectation_k_F <- 3 * n * (R - 1) / (3 * n * (R - 1) - 4)
  a <- 2 * (3 * R - 1) * expectation_k_F ** 2
  b <- 3 * (n - 2) * (R - 1) * (R ** 2) * (1 - 0.05 * R)^2
  phi <- (1+sqrt(5))/2
  if(M > 0){
    return(phi * ((1 + M) ** 4) * (a / b))
  }
  else{
    return(a/b)
  }

}

