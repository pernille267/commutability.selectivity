#' Approximate expectation of k given number of clinical samples and replicates
#'
#' @param n Integer - The number of unique clinical samples
#' @param R Integer - The number of replicated measurements done on each clinical sample
#'
#' @return A value that approximates the true expectation of k
#' @export
#'
#' @examples approximate_expectation_k(25,3)
approximate_expectation_k <- function(n, R){
  a <- 3 * n * (R - 1) / (R * (1 - 0.05 * R)) / (3 * n * (R-1) - 4)
  b <- sin(2) / (R * (1 - 0.05 * R))
  return(a - b + 1)
}

