#' Simulate k values based on simulation function
#'
#' @param m The number of simulated k values for each row of \code{parameters}
#' @param method The name of the function of type mc_* that you want to simulate k with. Allowed methods are \code{sdwinsp} (default), \code{sdwinsp2}, \code{sdwdnsp1}, \code{sdwdnsp2}, or 1, 2, 3, 4 respectively
#' @param parameters A data table or data frame where each row include a set of parameters. To generate a data table of unique combinations for a set of vectors of parameters it may be used to employ \code{CJ()} from the \code{data.table} package
#' @param to_dt If TRUE resulting k values are gathered into an data table
#'
#' @description A wrapper function for simulation of k values based on using either mc_* functions to simulate data. Note that all parameters not found in \code{parameters} and is relevant for the given \code{method} will be randomly generated or fixed. Special method dependent parameters such as xi, mmax and q will be fixed if not given in \code{parameters}
#'
#' @return A list or data table consisting of simulated k values based on generic parameters such as n, R, etc., and method dependent parameters. If \code{to_dt} is TRUE, the parameters found in \code{parameters} will be attached to the output data table
#' @export
#'
#' @examples simulate_k(m = 1e3, method = 3, parameters = data.frame(n = c(20,30), R = c(3, 3), p = c(0.05, 0.05), mmax = c(5,5)), to_dt = TRUE)

simulate_k <- function(m = 100, method = "sdwinsp", parameters, to_dt = TRUE){
  parameters_reference <- as.data.table(parameters)[,`parameter combination id`:=as.character(1:nrow(.SD))]
  parameters <- split(x = parameters_reference, f = 1:nrow(parameters_reference))
  if(any(method == "sdwinsp", method==1)){
    simulated_data <- lapply(X = parameters, FUN = mc_sdwinsp, m = m)
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
  else if(any(method == "sdwinsp2", method==2)){
    simulated_data <- lapply(X = parameters, FUN = mc_sdwinsp2, m = m)
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
  else if(any(method == "sdwdnsp1", method==3)){
    simulated_data <- lapply(X = parameters, FUN = mc_sdwdnsp1, m = m)
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
  else if(any(method == "sdwdnsp2", method==4)){
    simulated_data <- lapply(X = parameters, FUN = mc_sdwdnsp2, m = m)
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
  else{
    simulated_data <- lapply(X = parameters, FUN = mc_sdwinsp, m = m)
    simulated_data <- lapply(X = simulated_data, FUN = rbindlist, idcol = "simulation id internal")
    simulated_k_values <- lapply(X = simulated_data, FUN = estimate_k_over_groups, groups = "simulation id internal")
    if(to_dt){
      simulated_k_values <- rbindlist(l = simulated_k_values, idcol = "parameter combination id")
      message("input method is not recongnized. 'sdwinsp' is used. If you do not want this, you need to make sure your input matches the requirement. See ?simulate_k for details")
      return(merge.data.table(x = parameters_reference, y = simulated_k_values, by = "parameter combination id"))
    }
    else{
      message("input method is not recongnized. 'sdwinsp' is used. If you do not want this, you need to make sure your input matches the requirement. See ?simulate_k for details")
      return(list(k_values = simulated_k_values, parameters_lookup = parameters_reference))
    }
  }
}


