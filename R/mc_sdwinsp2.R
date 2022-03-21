#' Simulate several data sets with identical non-selectivity profiles using parallel computation
#'
#' @param parameter_row A data frame or data table with only 1 row. Should contain a subset of simulation parameters. Possible simulation parameters are n, R, cvx, cvy, ci_lwr and ci_upr, xi0, xi
#' @param m Integer - The number of data sets to be simulated
#' @param n Integer or 'r' - number of unique samples
#' @param R Integer or 'r' - Number of unique replicated measurements on each sample
#' @param cv_vec Vector with two elements or 'r' - MS CVs in decimal
#' @param ci Vector with two elements or 'r' - Concentration interval
#' @param xi0 A number signifying the relative starting value to the base. Default is 1 and should be smaller than \code{xi}
#' @param xi A number signifying the relative ending value relative to the base. Default is 2 and should be larger than \code{xi0}
#' @param parallel Logical - Should computations be performed in parallel or not. Default is \code{FALSE} and is not recommended unless m > 1e5
#' @param progress_bar Logical - Should a progress bar tracking the progress of the simulations be displayed? Default is \code{TRUE}
#'
#' @details Simulates m data sets with identical selectivity profiles, but with heteroscedasticity. The ratio of parameters \code{xi} and \code{xi0} indicates the relative increase of standard deviation of measurement errors from the lower end of the concentration range to the upper.
#'
#' @return A list containing m simulated data tables with 4 columns, two of which are measured results from two MSs in comparison. In addition, n * R rows.
#' @export
#'
#' @examples mc_sdwinsp2(m = 50, xi0 = 1, xi = 3)
mc_sdwinsp2 <- function(parameter_row = c(1), m = 1e2, n = "r", R = "r", cv_vec = "r", ci = "r", xi0 = 1, xi = 2, parallel = FALSE, progress_bar = TRUE){
  if(parallel){
    cl <- makeCluster(detectCores() - 1)
    clusterExport(cl=cl,varlist=c("sdwinsp2"))
    clusterEvalQ(cl=cl,expr=library(commutability.selectivity))
    if(progress_bar){
      mcs <- pblapply(cl = cl, X = as.list(1:m), FUN = function(x) sdwinsp2(n = n, R = R, cv_vec = cv_vec, ci = ci, xi0 = xi0, xi = xi, parameter_row = parameter_row))
    }
    else{
      mcs <- parLapply(cl = cl, X = as.list(1:m), fun = function(x) sdwinsp2(n = n, R = R, cv_vec = cv_vec, ci = ci, xi0 = xi0, xi = xi, parameter_row = parameter_row))
    }
    on.exit(expr = stopCluster(cl=cl), add = TRUE)
    return(mcs)
  }
  else{
    if(progress_bar){
      mcs <- pblapply(X = as.list(1:m), FUN = function(x) sdwinsp2(n = n, R = R, cv_vec = cv_vec, ci = ci, xi0 = xi0, xi = xi, parameter_row = parameter_row))
    }
    else{
      mcs <- lapply(X = as.list(1:m), FUN = function(x) sdwinsp2(n = n, R = R, cv_vec = cv_vec, ci = ci, xi0 = xi0, xi = xi, parameter_row = parameter_row))
    }
  }
}

