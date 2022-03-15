#' Simulate several data sets with identical non-selectivity profiles using parallel computation
#'
#' @param parameter_row A data frame or data table with only 1 row. Should contain a subset of simulation parameters. Possible simulation parameters are n, R, cvx, cvy, ci_lwr and ci_upr
#' @param m Integer - The number of data sets to be simulated
#' @param n Integer or 'r' - number of unique samples
#' @param R Integer or 'r' - Number of unique replicated measurements on each sample
#' @param cv_vec Vector with two elements or 'r' - MS CVs in decimal
#' @param ci Vector with two elements or 'r' - Concentration interval
#'
#' @details this function is recommended to use when simulating more than 1000 data sets, or else it will be faster to just use standard sapply
#'
#' @return A list containing m data tables with 4 columns, two of which are measured results from two MSs in comparison. In addition, n * R rows.
#' @export
#'
#' @examples mc_sdwinsp(m = 50)
mc_sdwinsp <- function(parameter_row = c(1), m = 1e2, n = "r", R = "r", cv_vec = "r", ci = "r"){
  cl <- makeCluster(detectCores() - 1)
  clusterExport(cl=cl,varlist=c("sdwinsp"))
  clusterEvalQ(cl=cl,expr=library(commutability.selectivity))
  mcs <- pbsapply(cl = cl, X = 1:m, FUN = sdwinsp, list(n = n, R = R, cv_vec = cv_vec, ci = ci, parameter_row = parameter_row), simplify = FALSE)
  on.exit(expr = stopCluster(cl=cl))
  return(mcs)
}



