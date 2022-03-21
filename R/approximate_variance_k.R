#' Approximate variance of k given number of clinical samples and replicates
#'
#' @param n Integer - The number of unique clinical samples
#' @param R Integer - The number of replicated measurements done on each clinical sample
#'
#' @return A value that approximates the true variance of k
#' @export
#'
#' @examples approximate_variance_k(25, 3)

approximate_variance_k <- function(n, R){
  expectation_k_F <- 3 * n * (R - 1) / (3 * n * (R - 1) - 4)
  a <- 2 * (3 * R - 1) * expectation_k_F ** 2
  b <- 3 * (n - 2) * (R - 1) * (R ** 2) * (1 - 0.05 * R)^2
  return(a/b)
}

