#' Simulation of data with identical ns profiles
#'
#' @param n Integer or 'r' - number of unique samples
#' @param R Integer or 'r' - number of unique samples
#' @param cv_vec vector with two elements or 'r' - MS CVs in decimal
#' @param ci vector with two elements or 'r' - Concentration interval
#' @param xi0 number that corresponds to starting SD respect to base SD
#' @param xi number that corresponds to ending SD respect to starting SD
#' @param parameter_row A data frame or data table with only 1 row. Should contain a subset of simulation parameters
#'
#' @description Simulation of data without non-selectivity but with heteroscedasticity. To be used in simulations
#'
#' @details Combine with lapply to simulate more than one data table. This function should only be used by advanced users. It is not user friendly for those not familiar with R
#'
#' @return A data table with 4 columns, two of which are measurement results from two MSs in comparison. In addition, n * R rows
#' @export
#'
#' @examples sdwinsp2()

sdwinsp2 <- function(n = "r", R = "r", cv_vec = "r", ci = "r", xi0 = 1, xi = 2, parameter_row = c(1)){
  cidomain <- sample(x=c(1e-1,1,2,3,5),size=1,prob=c(0.30,0.30,0.20,0.10,0.10))
  if(xi0==xi){xi = xi0 + 1e-5}
  n <- ifelse("n" %in% names(parameter_row), parameter_row$n, ifelse(n == "r", min(max(20, rpois(1, lambda = 25)), 30), n))
  R <- ifelse("R" %in% names(parameter_row), parameter_row$R, ifelse(R == "r", sample.int(5,1,prob=c(0.05,0.75,0.15,0.025,0.025))+1, R))
  cvx <- ifelse("cvx" %in% names(parameter_row), parameter_row$cvx, ifelse(cv_vec == "r", rbeta(1,2,5)/10, cv_vec[1]))
  cvy <- ifelse("cvy" %in% names(parameter_row), parameter_row$cvy, ifelse(cv_vec == "r", rbeta(1,2,5)/10, cv_vec[2]))
  ci_lwr <- ifelse("ci_lwr" %in% names(parameter_row), parameter_row$ci_lwr, ifelse(ci == "r", runif(1,0,1e2) * cidomain, ci[1]))
  ci_upr <- ifelse("ci_upr" %in% names(parameter_row), parameter_row$ci_upr, ifelse(ci == "r", ci_lwr + runif(1,0,1e2) * cidomain, ci[2]))
  tau <- runif(n = n, min = ci_lwr, max = ci_upr)
  tau <- rep(tau, each = R)
  xis <- rep(seq(from = min(xi0, xi), to = max(xi,xi0), length.out = n),each=R)
  rids <- rep(1:R, times = n)
  sids <- rep(1:n, each = R)
  sdx <- mean(c(ci_lwr,ci_upr)) * cvx
  sdy <- mean(c(ci_lwr,ci_upr)) * cvy
  dt <- setorder(x = data.table(SampleID=sids,ReplicateID=rids,tau=tau), tau)[,xis:=xis]
  dt <- data.table(SampleID=sids,ReplicateID=rids,tau=tau)[,.(ReplicateID=ReplicateID,MP_A=rnorm(R,tau,sdy),MP_B=rnorm(R,tau,sdx)),by=SampleID]
  return(dt)
}



