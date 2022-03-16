#' Simulate several data sets with different non-selectivity profiles defined by setting 1 using parallel computation
#'
#' @param parameter_row A data frame or data table with only 1 row. Should contain a subset of simulation parameters. Parameters that may be included are n, R, cvx, cvy, ci_lwr, ci_upr, mmax and p
#' @param m Integer - The number of data sets to be simulated
#' @param n Integer or 'r' - number of unique samples
#' @param R Integer or 'r' - Number of unique replicated measurements on each sample
#' @param cv_vec Vector with two elements or 'r' - MS CVs in decimal
#' @param ci Vector with two elements or 'r' - Concentration interval
#' @param mmax A number signifying the maximum multiple of relocation of clinical samples affected by differences in selectivity
#' @param p A number signifying the average proportion of clinical samples to be relocated
#' @param parallel allow parallelization
#'
#' @details this function is recommended to use when simulating more than 350 data sets, or else it will be faster to just use standard \code{sapply()}
#'
#' @return A list of simulated data tables with 4 columns, where two of which are measurement results from two MSs in comparison having differences in selectivity profiles. In addition, n * R rows
#' @export
#'
#' @examples mc_sdwdnsp1()
mc_sdwdnsp1 <- function(parameter_row = c(1), m = 1e2, n = "r", R = "r", cv_vec = "r", ci = "r", mmax = 2.5, p = 0.05, parallel = TRUE){
  if(parallel){
    cl <- makeCluster(detectCores() - 1)
    clusterExport(cl = cl, varlist = c("sdwdnsp1"))
    clusterEvalQ(cl = cl, expr = library(commutability.selectivity))
    mcs <- pbsapply(cl = cl, X = 1:m, FUN = sdwdnsp1, list(n = n, R = R, cv_vec = cv_vec, ci = ci, mmax = mmax, p = p, parameter_row = parameter_row), simplify = FALSE)
    on.exit(expr = stopCluster(cl=cl), add = TRUE)
    return(mcs)
  }
  else{
    mcs <- sapply(X = 1:m, FUN = sdwdnsp1, list(n = n, R = R, cv_vec = cv_vec, ci = ci, mmax = mmax, p = p, parameter_row = parameter_row), simplify = FALSE)
    return(mcs)
  }


}
