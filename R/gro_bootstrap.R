#' bootstrap estimates
#'
#' @param data some data
#' @param filter_var some value referring to missing values
#' @param B how many bootstrap replicates should be used
#' @param alpha significance level
#' @param statistic Which statistic should we do inference on. A function that exists in R taking a numeric vector. median is default
#' @param calculate What do we want to calculate for the statistic considered. EIther \code{confidence interval}, \code{variance}, \code{standard deviation} or \code{coefficient of variation}
#'
#' @return values for each considered group
#' @export
#'
#' @examples print(1)
gro_bootstrap <- function(data, filter_var = "NO", B = 1000, alpha = 0.05, statistic = median, calculate = c("confidence interval","variance","standard deviation","coefficient of variation")){
  data <- setDT(data)[,lapply(X = .SD, FUN = stri_replace_all, replace = NA, fixed = filter_var)]
  data <- melt.data.table(data = setDT(data),
                          measure.vars = names(data),
                          variable.name = "MS",
                          value.name = "RD",
                          na.rm = TRUE)[,.(RD = as.numeric(RD)), by = MS]
  new_data <- split(x = data, by = "MS")
  cl <- makeCluster(detectCores() - 1)
  clusterEvalQ(cl = cl, expr = library(commutability.selectivity))
  new_data <- pblapply(X = new_data,
                       FUN = function(x) gro_resample(data = x, B = B),
                       cl = cl)
  new_data <- pblapply(X = new_data,
                       FUN = function(x) x[,lapply(X=.SD,FUN=statistic),.SDcols="RD",by="RSTEP"])
  if(calculate[1]=="confidence interval"){
    calculation <- pblapply(X = new_data,
                            FUN = function(x) x[,.(center = mean(RD),
                                                   lwr = quantile(x = RD, probs = alpha/2),
                                                   upr = quantile(x = RD, probs = 1 - alpha/2))])
  }
  else if(calculate[1]=="variance"){
    calculation <- lapply(X = new_data,
                            FUN = function(x) x[,.(center = mean(RD), variance = var(RD))])
  }

  else if(calculate[1]=="standard deviation"){
    calculation <- lapply(X = new_data,
                            FUN = function(x) x[,.(center = mean(RD), `standard deviation` = var(RD))])
  }

  else if(calculate[1]=="coefficient of variation"){
    calculation <- lapply(X = new_data,
                            FUN = function(x) x[,.(center = mean(RD),
                                                   cv = sd(RD)/mean(RD),
                                                   robcv = mad(RD)/median(RD))])
  }
  else{
    calculation <- lapply(X = new_data,
                            FUN = function(x) x[,.(center = mean(RD),
                                                   lwr = quantile(x = RD, probs = alpha/2),
                                                   upr = quantile(x = RD, probs = 1 - alpha/2))])
    message("no valid value for 'calculate' found. Accordingly, confidence interval is calculated for statistic instead")
  }
  calculation <- rbindlist(l = calculation, idcol = "Measurement system")
  on.exit(expr = stopCluster(cl = cl), add = TRUE)
  return(calculation)
}



