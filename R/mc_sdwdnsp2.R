#' Simulate several data sets with different non-selectivity profiles defined by setting 2 using parallel computation
#'
#' @param parameter_row A data frame or data table with only 1 row. Should contain a subset of simulation parameters. Parameters that may be included are n, R, cvx, cvy, ci_lwr, ci_upr, mmax, qmin and  qmax
#' @param m Integer - The number of data sets to be simulated
#' @param n Integer or 'r' - number of unique samples
#' @param R Integer or 'r' - Number of unique replicated measurements on each sample
#' @param cv_vec Vector with two elements or 'r' - MS CVs in decimal
#' @param ci Vector with two elements or 'r' - Concentration interval
#' @param mmax A number signifying the maximum multiple of relocation of clinical samples affected by differences in selectivity
#' @param q A vector with two elements signifying the lower and upper quantile boundaries where clinical samples are relocated
#'
#' @details this function is recommended to use when simulating more than 350 data sets, or else it will be faster to just use standard \code{sapply()}. This function is not available for other users than Windows.
#'
#' @return A list of simulated data tables with 4 columns, where two of which are measurement results from two MSs in comparison having differences in selectivity profiles defined by setting 2. In addition, n * R rows
#' @export
#'
#' @examples mc_sdwdnsp2(m = 50)
mc_sdwdnsp2 <- function(parameter_row = c(1), m = 1e2, n = "r", R = "r", cv_vec = "r", ci = "r", mmax = 2.5, q = c(0,0.25)){
  cl <- makeCluster(detectCores() - 1)
  clusterExport(cl = cl, varlist = c("sdwdnsp2"))
  clusterEvalQ(cl = cl, expr = library(commutability.selectivity))
  mcs <- pbsapply(cl = cl, X = 1:m, FUN = sdwdnsp2, list(n = n, R = R, cv_vec = cv_vec, ci = ci, mmax = mmax, q = q, parameter_row = parameter_row), simplify = FALSE)
  on.exit(expr = stopCluster(cl=cl))
  return(mcs)
}
