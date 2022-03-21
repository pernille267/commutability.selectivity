#' Simulate k values based on simulation function
#'
#' @param m The number of simulated k values for each row of \code{parameters}
#' @param method The name of the function of type mc_* that you want to simulate k with. Allowed methods are \code{sdwinsp} (default), \code{sdwinsp2}, \code{sdwdnsp1}, \code{sdwdnsp2}, \code{sdwdnsp3}, or 1, 2, 3, 4, 5 respectively
#' @param parameters A data table or data frame where each row include a set of parameters. To generate a data table of unique combinations for a set of vectors of parameters it may be used to employ \code{CJ()} from the \code{data.table} package
#' @param to_dt If TRUE resulting k values are gathered into an data table
#' @param parallel If TRUE, parallel computing will be used for each parameter combination instead of each iterations per parameter combination
#'
#' @description A wrapper function for simulation of k values based on using either mc_* functions to simulate data. Note that all parameters not found in \code{parameters} and is relevant for the given \code{method} will be randomly generated or fixed. Special method dependent parameters such as xi, mmax and q will be fixed if not given in \code{parameters}
#'
#' @return A list or data table consisting of simulated k values based on generic parameters such as n, R, etc., and method dependent parameters. If \code{to_dt} is TRUE, the parameters found in \code{parameters} will be attached to the output data table
#' @export
#'
#' @examples simulate_k(m = 1e3, method = 3, parameters = data.frame(n = c(20,30), R = c(3, 3), p = c(0.05, 0.05), mmax = c(5,5)), to_dt = TRUE)

simulate_k <- function(m = 100, method = "sdwinsp", parameters, to_dt = TRUE, parallel = TRUE){
  parameters_reference <- as.data.table(parameters)[,`parameter combination id`:=as.character(1:nrow(.SD))]
  parameters <- split(x = parameters_reference, f = 1:nrow(parameters_reference))
  if(any(method == "sdwinsp", method==1)){
    if(parallel){
      cl <- makeCluster(detectCores() - 1)
      clusterEvalQ(cl = cl, expr = library(commutability.selectivity))
      message("simulating data ...")
      simulated_data <- pblapply(X = parameters, FUN = mc_sdwinsp, m = m, cl = cl, parallel = all(!parallel,length(parameters)<10,m>5e4), progress_bar = FALSE)
      simulated_data <- lapply(X = simulated_data, FUN = rbindlist, idcol = "simulation id internal")
      message("calculating k values based on simulated data ...")
      simulated_k_values <- pblapply(X = simulated_data, FUN = estimate_k_over_groups, groups = "simulation id internal", cl = cl)
      if(to_dt){
        simulated_k_values <- rbindlist(l = simulated_k_values, idcol = "parameter combination id")
        on.exit(expr = stopCluster(cl = cl), add = TRUE)
        return(merge.data.table(x = parameters_reference, y = simulated_k_values, by = "parameter combination id"))
      }
      else{
        on.exit(expr = stopCluster(cl = cl), add = TRUE)
        return(list(k_values = simulated_k_values, parameters_lookup = parameters_reference))
      }
    }
    else{
      simulated_data <- lapply(X = parameters, FUN = mc_sdwinsp, m = m, parallel = all(!parallel,length(parameters)<10,m>5e4), progress_bar = FALSE)
      simulated_data <- lapply(X = simulated_data, FUN = rbindlist, idcol = "simulation id internal")
      simulated_k_values <- lapply(X = simulated_data, FUN = estimate_k_over_groups, groups = "simulation id internal")
      if(to_dt){
        simulated_k_values <- rbindlist(l = simulated_k_values, idcol = "parameter combination id")
        return(merge.data.table(x = parameters_reference, y = simulated_k_values, by = "parameter combination id"))
      }
      else{
        return(list(k_values = simulated_k_values, parameters_lookup = parameters_reference))
      }
    }
  }
  else if(any(method == "sdwinsp2", method==2)){
    if(parallel){
      cl <- makeCluster(detectCores() - 1)
      clusterEvalQ(cl = cl, expr = library(commutability.selectivity))
      message("simulating data ...")
      simulated_data <- pblapply(X = parameters, FUN = mc_sdwinsp2, m = m, cl = cl, parallel = all(!parallel,length(parameters)<10,m>5e4), progress_bar = FALSE)
      simulated_data <- lapply(X = simulated_data, FUN = rbindlist, idcol = "simulation id internal")
      message("calculating k values based on simulated data ...")
      simulated_k_values <- pblapply(X = simulated_data, FUN = estimate_k_over_groups, groups = "simulation id internal", cl = cl)
      if(to_dt){
        simulated_k_values <- rbindlist(l = simulated_k_values, idcol = "parameter combination id")
        on.exit(expr = stopCluster(cl = cl), add = TRUE)
        return(merge.data.table(x = parameters_reference, y = simulated_k_values, by = "parameter combination id"))
      }
      else{
        on.exit(expr = stopCluster(cl = cl), add = TRUE)
        return(list(k_values = simulated_k_values, parameters_lookup = parameters_reference))
      }
    }
    else{
      simulated_data <- lapply(X = parameters, FUN = mc_sdwinsp2, m = m, parallel = all(!parallel,length(parameters)<10,m>5e4), progress_bar = FALSE)
      simulated_data <- lapply(X = simulated_data, FUN = rbindlist, idcol = "simulation id internal")
      simulated_k_values <- lapply(X = simulated_data, FUN = estimate_k_over_groups, groups = "simulation id internal")
      if(to_dt){
        simulated_k_values <- rbindlist(l = simulated_k_values, idcol = "parameter combination id")
        return(merge.data.table(x = parameters_reference, y = simulated_k_values, by = "parameter combination id"))
      }
      else{
        return(list(k_values = simulated_k_values, parameters_lookup = parameters_reference))
      }
    }
  }
  else if(any(method == "sdwdnsp1", method==3)){
    if(parallel){
      cl <- makeCluster(detectCores() - 1)
      clusterEvalQ(cl = cl, expr = library(commutability.selectivity))
      message("simulating data with differences in selectivity defined by setting 1 ...")
      simulated_data <- pblapply(X = parameters, FUN = mc_sdwdnsp1, m = m, cl = cl, parallel = all(!parallel,length(parameters)<10,m>5e4), progress_bar = FALSE)
      simulated_data <- lapply(X = simulated_data, FUN = rbindlist, idcol = "simulation id internal")
      message("calculating k values based on simulated data ...")
      simulated_k_values <- pblapply(X = simulated_data, FUN = estimate_k_over_groups, groups = "simulation id internal", cl = cl)
      if(to_dt){
        simulated_k_values <- rbindlist(l = simulated_k_values, idcol = "parameter combination id")
        on.exit(expr = stopCluster(cl = cl), add = TRUE)
        return(merge.data.table(x = parameters_reference, y = simulated_k_values, by = "parameter combination id"))
      }
      else{
        on.exit(expr = stopCluster(cl = cl), add = TRUE)
        return(list(k_values = simulated_k_values, parameters_lookup = parameters_reference))
      }
    }
    else{
      simulated_data <- lapply(X = parameters, FUN = mc_sdwdnsp1, m = m, parallel = all(!parallel,length(parameters)<10,m>5e4), progress_bar = FALSE)
      simulated_data <- lapply(X = simulated_data, FUN = rbindlist, idcol = "simulation id internal")
      simulated_k_values <- lapply(X = simulated_data, FUN = estimate_k_over_groups, groups = "simulation id internal")
      if(to_dt){
        simulated_k_values <- rbindlist(l = simulated_k_values, idcol = "parameter combination id")
        return(merge.data.table(x = parameters_reference, y = simulated_k_values, by = "parameter combination id"))
      }
      else{
        return(list(k_values = simulated_k_values, parameters_lookup = parameters_reference))
      }
    }
  }
  else if(any(method == "sdwdnsp2", method==4)){
    if(parallel){
      cl <- makeCluster(detectCores() - 1)
      clusterEvalQ(cl = cl, expr = library(commutability.selectivity))
      message("simulating data with differences in selectivity defined by setting 2 ...")
      simulated_data <- pblapply(X = parameters, FUN = mc_sdwdnsp2, m = m, cl = cl, parallel = all(!parallel,length(parameters)<10,m>5e4), progress_bar = FALSE)
      simulated_data <- lapply(X = simulated_data, FUN = rbindlist, idcol = "simulation id internal")
      message("calculating k values based on simulated data ...")
      simulated_k_values <- pblapply(X = simulated_data, FUN = estimate_k_over_groups, groups = "simulation id internal", cl = cl)
      if(to_dt){
        simulated_k_values <- rbindlist(l = simulated_k_values, idcol = "parameter combination id")
        on.exit(expr = stopCluster(cl = cl), add = TRUE)
        return(merge.data.table(x = parameters_reference, y = simulated_k_values, by = "parameter combination id"))
      }
      else{
        on.exit(expr = stopCluster(cl = cl), add = TRUE)
        return(list(k_values = simulated_k_values, parameters_lookup = parameters_reference))
      }
    }
    else{
      simulated_data <- lapply(X = parameters, FUN = mc_sdwdnsp2, m = m, parallel = all(!parallel,length(parameters)<10,m>5e4), progress_bar = FALSE)
      simulated_data <- lapply(X = simulated_data, FUN = rbindlist, idcol = "simulation id internal")
      simulated_k_values <- lapply(X = simulated_data, FUN = estimate_k_over_groups, groups = "simulation id internal")
      if(to_dt){
        simulated_k_values <- rbindlist(l = simulated_k_values, idcol = "parameter combination id")
        return(merge.data.table(x = parameters_reference, y = simulated_k_values, by = "parameter combination id"))
      }
      else{
        return(list(k_values = simulated_k_values, parameters_lookup = parameters_reference))
      }
    }
  }
  else if(any(method == "sdwdnsp3", method==5)){
    if(parallel){
      cl <- makeCluster(detectCores() - 1)
      clusterEvalQ(cl = cl, expr = library(commutability.selectivity))
      message("simulating data with differences in selectivity defined by setting 3 ...")
      simulated_data <- pblapply(X = parameters, FUN = mc_sdwinsp, m = m, cl = cl, parallel = all(!parallel,length(parameters)<10,m>5e4), progress_bar = FALSE)
      simulated_data <- lapply(X = simulated_data, FUN = rbindlist, idcol = "simulation id internal")
      message("calculating k values based on simulated data ...")
      simulated_k_values <- pblapply(X = simulated_data, FUN = estimate_k_over_groups, groups = "simulation id internal", cl = cl)
      if(to_dt){
        simulated_k_values <- rbindlist(l = simulated_k_values, idcol = "parameter combination id")
        simulated_k_values <- simulated_k_values
        on.exit(expr = stopCluster(cl = cl), add = TRUE)
        return(merge.data.table(x = parameters_reference, y = simulated_k_values, by = "parameter combination id")[,.(`parameter combination id`, n = n, R = R, M = M, k = k*(1+M)**2)])
      }
      else{
        on.exit(expr = stopCluster(cl = cl), add = TRUE)
        return(list(k_values = simulated_k_values, parameters_lookup = parameters_reference))
      }
    }
    else{
      simulated_data <- lapply(X = parameters, FUN = mc_sdwinsp, m = m, parallel = all(!parallel,length(parameters)<10,m>5e4), progress_bar = FALSE)
      simulated_data <- lapply(X = simulated_data, FUN = rbindlist, idcol = "simulation id internal")
      simulated_k_values <- lapply(X = simulated_data, FUN = estimate_k_over_groups, groups = "simulation id internal")
      if(to_dt){
        simulated_k_values <- rbindlist(l = simulated_k_values, idcol = "parameter combination id")
        simulated_k_values <- simulated_k_values
        return(merge.data.table(x = parameters_reference, y = simulated_k_values, by = "parameter combination id")[,.(`parameter combination id`, n = n, R = R, M = M, k = k*(1+M)**2)])
      }
      else{
        return(list(k_values = simulated_k_values, parameters_lookup = parameters_reference))
      }
    }
  }
  else{
    if(parallel){
      cl <- makeCluster(detectCores() - 1)
      clusterEvalQ(cl = cl, expr = library(commutability.selectivity))
      message("simulating data from sdwinsp because no valid method are chosen ...")
      simulated_data <- pblapply(X = parameters, FUN = mc_sdwinsp, m = m, cl = cl, parallel = all(!parallel,length(parameters)<10,m>5e4), progress_bar = FALSE)
      simulated_data <- lapply(X = simulated_data, FUN = rbindlist, idcol = "simulation id internal")
      message("calculating k values based on simulated data ...")
      simulated_k_values <- pblapply(X = simulated_data, FUN = estimate_k_over_groups, groups = "simulation id internal", cl = cl)
      if(to_dt){
        simulated_k_values <- rbindlist(l = simulated_k_values, idcol = "parameter combination id")
        on.exit(expr = stopCluster(cl = cl), add = TRUE)
        return(merge.data.table(x = parameters_reference, y = simulated_k_values, by = "parameter combination id"))
      }
      else{
        on.exit(expr = stopCluster(cl = cl), add = TRUE)
        return(list(k_values = simulated_k_values, parameters_lookup = parameters_reference))
      }
    }
    else{
      simulated_data <- lapply(X = parameters, FUN = mc_sdwinsp, m = m, parallel = all(!parallel,length(parameters)<10,m>5e4), progress_bar = FALSE)
      simulated_data <- lapply(X = simulated_data, FUN = rbindlist, idcol = "simulation id internal")
      simulated_k_values <- lapply(X = simulated_data, FUN = estimate_k_over_groups, groups = "simulation id internal")
      if(to_dt){
        simulated_k_values <- rbindlist(l = simulated_k_values, idcol = "parameter combination id")
        return(merge.data.table(x = parameters_reference, y = simulated_k_values, by = "parameter combination id"))
      }
      else{
        return(list(k_values = simulated_k_values, parameters_lookup = parameters_reference))
      }
    }
  }
}


