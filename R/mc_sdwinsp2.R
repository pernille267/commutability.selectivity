#' Simulate several data sets with identical non-selectivity profiles using parallel computation
#'
#' @param parameter_row A data frame or data table with only 1 row. Should contain a subset of simulation parameters. Possible simulation parameters are n, R, cvx, cvy, ci_lwr and ci_upr, xi0, xi
#' @param m Integer - The number of data sets to be simulated
#' @param n Integer or 'r' - number of unique samples
#' @param R Integer or 'r' - Number of unique replicated measurements on each sample
#' @param cv_vec Vector with two elements or 'r' - MS CVs in decimal
#' @param ci Vector with two elements or 'r' - Concentration interval
#' @param xi0 A number signifying the relative starting value to the base. Default is 1
#' @param xi A number signifying the relative ending value relative to the base. Default is 2
#'
#' @details this function is recommended to use when simulating more than 500 data sets, or else it will be faster to just use standard sapply
#'
#' @return A list containing m simulated data tables with 4 columns, two of which are measured results from two MSs in comparison. In addition, n * R rows.
#' @export
#'
#' @examples mc_sdwinsp2(m = 50, xi0 = 1, xi = 3)
mc_sdwinsp2 <- function(parameter_row = c(1), m = 1e2, n = "r", R = "r", cv_vec = "r", ci = "r", xi0 = 1, xi = 2){
  cl <- makeCluster(detectCores() - 1)
  clusterExport(cl=cl,varlist=c("sdwinsp2"))
  clusterEvalQ(cl=cl,expr=library(commutability.selectivity))
  mcs <- pbsapply(cl = cl, X = 1:m, FUN = sdwinsp2, list(n = n, R = R, cv_vec = cv_vec, ci = ci, xi0 = xi0, xi = xi, parameter_row = parameter_row), simplify = FALSE)
  on.exit(expr = stopCluster(cl=cl))
  return(mcs)
}
