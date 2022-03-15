#' Simulation of data with identical ns profiles
#'
#' @param n Integer or 'r' - number of unique samples
#' @param R Integer or 'r' - number of unique samples
#' @param cv_vec vector with two elements or 'r' - MS CVs in decimal
#' @param ci vector with two elements or 'r' - Concentration interval
#' @param parameter_row A data frame or data table with only 1 row. Should contain a subset of simulation parameters
#'
#' @description Simulation of data with identical selectivity selectivity profiles. To be used in simulations. Setting either simulation parameters to 'r' will invoke random generation of these numbers. Note that specifying parameters both outside and inside of parameter_row will make those parameters outside being overridden.
#'
#' @details Combine with \code{lapply()} or \code{sapply()} to simulate more than one data table. For parallel computing use \code{mc_sdwinsp()}. This function should only be used by advanced users. It is not user friendly for those not familiar with R
#'
#' @return A data table with 4 columns, two of which are measurement results from two MSs in comparison. In addition, n * R rows
#' @export
#'
#' @examples sdwinsp()
sdwinsp <- function(n = "r", R = "r", cv_vec = "r", ci = "r", parameter_row = c(1)){
  cidomain <- sample(x=c(1e-1,1,2,3,5),size=1,prob=c(0.30,0.30,0.20,0.10,0.10))

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
  else if(ci == "r"){
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
  else if(ci == "r"){
    ci_upr <- ci_lwr + runif(1,0,1e2) * cidomain
  }
  else if(is.double(ci[2])){
    ci_upr <- ci[2]
  }
  else{
    ci_upr <- 100
  }

  tau <- runif(n = n, min = ci_lwr, max = ci_upr)
  tau <- rep(tau, each = R)
  rids <- rep(1:R, times = n)
  sids <- rep(1:n, each = R)
  sdx <- sum(c(ci_lwr,ci_upr))/2 * cvx
  sdy <- sum(c(ci_lwr,ci_upr))/2 * cvy
  dt <- list("SampleID"=sids, "ReplicateID"=rids, "tau"=tau)
  dt <- setDT(x = dt)[,.(ReplicateID=ReplicateID,MP_A=rnorm(R,tau,sdy),MP_B=rnorm(R,tau,sdx)),by=SampleID]
  return(dt)
}



