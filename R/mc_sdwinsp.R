#' Simulate several data sets with identical non-selectivity profiles using parallel computation
#'
#' @param parameter_row A data frame or data table with only 1 row. Should contain a subset of simulation parameters. Possible simulation parameters are n, R, cvx, cvy, ci_lwr and ci_upr
#' @param m Integer - The number of data sets to be simulated
#' @param n Integer or 'r' - number of unique samples
#' @param R Integer or 'r' - Number of unique replicated measurements on each sample
#' @param cv_vec Vector with two elements or 'r' - MS CVs in decimal
#' @param ci Vector with two elements or 'r' - Concentration interval
#' @param parallel Should we allow parallel computing to be performed when simulating from \code{sdwinsp()}. This is not recommended if the number of parameter combinations is larger than 10. Default value is FALSE
#' @param progress_bar Should a progress bar monitoring the simulation progress be displayed. Default is TRUE
#'
#' @details this function is recommended to use when simulating more than 1000 data sets, or else it will be faster to just use standard sapply
#'
#' @return A list containing m data tables with 4 columns, two of which are measured results from two MSs in comparison. In addition, n * R rows.
#' @export
#'
#' @examples mc_sdwinsp(m = 50)
mc_sdwinsp <- function(parameter_row = c(2), m = 1e2, n = "r", R = "r", cv_vec = "r", ci = "r", parallel = FALSE, progress_bar = TRUE){
  if(parallel){
    cl <- makeCluster(detectCores() - 1)
    clusterExport(cl = cl, varlist = c("sdwinsp"))
    clusterEvalQ(cl = cl, expr = library(commutability.selectivity))
    if(progress_bar){
      mcs <- pblapply(X = as.list(1:m), FUN = function(x) sdwinsp(n = n, R = R, cv_vec = cv_vec, ci = ci, parameter_row = parameter_row), cl = cl)
    }
    else{
      mcs <- parLapply(cl = cl, X = as.list(1:m), fun = function(x) sdwinsp(n = n, R = R, cv_vec = cv_vec, ci = ci, parameter_row = parameter_row))
    }
    on.exit(expr = stopCluster(cl = cl), add = TRUE)
    return(mcs)
  }
  else{
    if(progress_bar){
      mcs <- pblapply(X = as.list(1:m), FUN = function(x) sdwinsp(n = n, R = R, cv_vec = cv_vec, ci = ci, parameter_row = parameter_row))
      return(mcs)
    }
    else{
      mcs <- lapply(X = as.list(1:m), FUN = function(x) sdwinsp(n = n, R = R, cv_vec = cv_vec, ci = ci, parameter_row = parameter_row))
      return(mcs)
    }
  }
}



