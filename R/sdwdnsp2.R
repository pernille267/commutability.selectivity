#' Simulation of homoscedastic data with different non-selectivity profiles defined by the second setting
#'
#' @param n Integer or 'r' - number of unique samples
#' @param R Integer or 'r' - number of unique replicated measurements of each sample
#' @param cv_vec vector with two elements or 'r' - MS CVs in decimal
#' @param ci vector with two elements or 'r' - Concentration interval
#' @param mmax A number signifying the maximum multiple of relocation of clinical samples affected by differences in selectivity
#' @param q A vector with two elements signifying the lower and upper quantile boundaries where clinical samples are relocated
#' @param parameter_row A data frame or data table with only 1 row. Should contain a subset of simulation parameters. Parameters that may be included are n, R, cvx, cvy, ci_lwr, ci_upr, mmax, qmin and qmax
#' @param include_relocated Should an additional column be attached that signify which CSs that were relocated? Default is FALSE
#'
#' @description Simulation of data with different selectivity profiles defined by second setting. To be used in simulations. Setting either simulation parameters to 'r' will invoke random generation of these numbers. Note that specifying parameters both outside and inside of parameter_row will make those parameters outside being overridden.
#'
#' @details Combine with \code{lapply()} or \code{sapply()} to simulate more than one data table. For parallel computing use \code{mc_sdwdnsp2()}. This function should only be used by advanced users. It is not user friendly for those not familiar with R
#'
#' @return A data table with 4 columns, two of which are measurement results from two MSs in comparison having differences in selectivity profiles defined by second setting. In addition, n * R rows
#' @export
#'
#' @examples sdwdnsp2()
sdwdnsp2 <- function(n = "r", R = "r", cv_vec = "r", ci = "r", mmax = 2.5, q = c(0,0.25), parameter_row = c(1), include_relocated = FALSE){

  cidomain <- sample(x = c(1e-1, 1, 2, 3, 5), size = 1, prob = c(0.30,0.30,0.20,0.10,0.10))

  if(any(names(parameter_row) == "n")){
    n <- parameter_row$n
  }
  else if(n == "r"){
    n <- min(max(20, rpois(1, lambda = 25)), 30)
  }
  else if(is.integer(round(n))){
    n <- n
  }
  else{
    n <- 25
  }

  if(any(names(parameter_row) == "R")){
    R <- parameter_row$R
  }
  else if(R == "r"){
    R <- sample.int(5,1,prob=c(0.05,0.75,0.15,0.025,0.025))+1
  }
  else if(is.integer(R)){
    R <- R
  }
  else{
    R <- 3
  }

  if(any(names(parameter_row) == "cvx")){
    cvx <- parameter_row$cvx
  }
  else if(cv_vec[1] == "r"){
    cvx <- rbeta(1,2,5)/10
  }
  else if(is.double(cv_vec[1])){
    cvx <- cv_vec[1]
  }
  else{
    cvx <- 0.05
  }

  if(any(names(parameter_row) == "cvy")){
    cvy <- parameter_row$cvy
  }
  else if(cv_vec[1] == "r"){
    cvy <- rbeta(1,2,5)/10
  }
  else if(is.double(cv_vec[2])){
    cvy <- cv_vec[2]
  }
  else{
    cvy <- 0.05
  }

  if(any(names(parameter_row) == "ci_lwr")){
    ci_lwr <- parameter_row$ci_lwr
  }
  else if(ci[1] == "r"){
    ci_lwr <- runif(1,0,1e2) * cidomain
  }
  else if(is.double(ci[1])){
    ci_lwr <- ci[1]
  }
  else{
    ci_lwr <- 100
  }

  if(any(names(parameter_row) == "ci_upr")){
    ci_upr <- parameter_row$ci_upr
  }
  else if(ci[1] == "r"){
    ci_upr <- ci_lwr + runif(1,0,1e2) * cidomain
  }
  else if(is.double(ci[2])){
    ci_upr <- ci[2]
  }
  else{
    ci_upr <- 100
  }

  if(any(names(parameter_row) == "mmax")){
    mmax <- parameter_row$mmax
  }
  else if(is.double(mmax)){
    mmax <- mmax
  }
  else{
    mmax <- 2.5
  }

  if(any(names(parameter_row) == "qmin")){
    qmin <- parameter_row$qmin
  }
  else if(all(is.double(q))){
    q <- q
  }
  else{
    q <- c(0,0.25)
  }

  if(any(names(parameter_row) == "qmax")){
    qmax <- parameter_row$qmax
    q <- c(qmin,qmax)
  }
  else if(all(is.double(q))){
    q <- q
  }
  else{
    q <- c(0,0.25)
  }

  tau <- runif(n = n, min = ci_lwr, max = ci_upr)
  tau <- rep(tau, each = R)
  rids <- rep(1:R, times = n)
  sids <- rep(1:n, each = R)
  sdx <- (sum(c(ci_lwr,ci_upr))/2) * cvx
  sdy <- (sum(c(ci_lwr,ci_upr))/2) * cvy
  dir <- sample(x=c(1,-1),size=1)
  dt <- list("SampleID" = sids, "ReplicateID" = rids, "tau" = tau)
  dt <- setDT(dt)[, influenced := tau <= quantile(x = tau, probs = q[2], na.rm = TRUE, names = FALSE) & tau >= quantile(x = tau, probs = q[1], na.rm = TRUE, names = FALSE)]
  dt <- dt[,.(ReplicateID = ReplicateID,
              influenced = influenced,
              Relocation = influenced * dir * runif(1, 0, mmax * sqrt(sdy**2 + sdx**2)) / sqrt(2),
              MP_A = rnorm(R,tau,sdy),
              MP_B = rnorm(R,tau,sdx)), by = SampleID]
  if(include_relocated){
    dt <- dt[,.(ReplicateID = ReplicateID,
                Relocated = influenced == TRUE,
                MP_A = MP_A + Relocation,
                MP_B = MP_B - Relocation), by = SampleID]
  }
  else{
    dt <- dt[,.(ReplicateID = ReplicateID,
                MP_A = MP_A + Relocation * influenced,
                MP_B = MP_B - Relocation * influenced), by = SampleID]
  }
  return(dt)
}


