#' XXX
#'
#' @param data data from GRO
#' @param filter_var which value represents missing value in data
#' @param alpha significance level for two-sided asymptotic z-interval
#'
#' @return returns something useful
#' @export
#'
#' @examples print(1)
gro_standard <- function(data, filter_var = "NO", alpha = 0.05){
  data <- setDT(data)[,lapply(X = .SD, FUN = stri_replace_all, replace = NA, fixed = filter_var)]
  data <- melt.data.table(data = setDT(data),
                          measure.vars = names(data),
                          variable.name = "MS",
                          value.name = "RD",
                          na.rm = TRUE)[,.(RD = as.numeric(RD)), by = MS]
  data <- split(x = data, by = "MS")
  normality <- lapply(X = data, FUN = function(x) data.table(normality_test_p.value = shapiro.test(x$RD)$p.value))
  normality <- rbindlist(l = normality, idcol = "Measurement system")
  CI_norm <- lapply(X = data, FUN = function(x) x[,width := qnorm(p = 1-alpha/2) * sd(RD)/sqrt(length(RD))])
  CI_norm <- lapply(X = CI_norm, FUN = function(x) unique(x[,.(median = median(RD), lwr = median(RD) - width, upr = median(RD) + width)]))
  CI_norm <- rbindlist(l = CI_norm, idcol = "Measurement system")
  all_standard_results <- merge.data.table(x = CI_norm, y = normality, by = "Measurement system")
  return(all_standard_results)
}




