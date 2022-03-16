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
#' @param parallel should we allow parallelization to be performed
#' @param progress_bar should a progress bar be displayed. Default is \code{TRUE}
#'
#' @description Performance; 100 - TT 2 sec TF 2 sec FT 0.25 sec FF 0.15 sec, 1,000 - TT 3.9 sec TF 2.3 sec FT 2.1 sec FF 1.7 sec, 10,000 TT 8.1 sec TF 6 sec FT 18 sec FF 17.5 sec
#' @details \code{parallel} should only be TRUE used if parameter_row has few rows, or else it will be faster to just use standard \code{sapply()} by setting \code{parallel} to FALSE. This function is not available for other users than Windows.
#'
#' @return A list of simulated data tables with 4 columns, where two of which are measurement results from two MSs in comparison having differences in selectivity profiles defined by setting 2. In addition, n * R rows
#' @export
#'
#' @examples mc_sdwdnsp2(m = 50)
mc_sdwdnsp2 <- function(parameter_row = c(1), m = 1e2, n = "r", R = "r", cv_vec = "r", ci = "r", mmax = 2.5, q = c(0,0.25), parallel = TRUE, progress_bar = TRUE){
  if(parallel){
    cl <- makeCluster(detectCores() - 1)
    clusterEvalQ(cl = cl, expr = library(commutability.selectivity))
    if(progress_bar){
      mcs <- pbsapply(cl = cl, X = 1:m, FUN = sdwdnsp2, n = n, R = R, cv_vec = cv_vec, ci = ci, mmax = mmax, q = q, parameter_row = parameter_row, simplify = FALSE)
    }
    else{
      mcs <- parSapply(cl = cl, X = 1:m, FUN = sdwdnsp2, n = n, R = R, cv_vec = cv_vec, ci = ci, mmax = mmax, q = q, parameter_row = parameter_row, simplify = FALSE)
    }
    on.exit(expr = stopCluster(cl=cl), add = TRUE)
    return(mcs)
  }
  else{
    if(progress_bar){
      mcs <- pbsapply(X = 1:m, FUN = sdwdnsp2, n = n, R = R, cv_vec = cv_vec, ci = ci, mmax = mmax, q = q, parameter_row = parameter_row, simplify = FALSE)
    }
    else{
      mcs <- sapply(X = 1:m, FUN = sdwdnsp2, n = n, R = R, cv_vec = cv_vec, ci = ci, mmax = mmax, q = q, parameter_row = parameter_row, simplify = FALSE)
    }
    return(mcs)
  }
}
