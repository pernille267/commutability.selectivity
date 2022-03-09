#' Estimate precision of measurement systems in comparison based on replicated measurements for LFDT data
#'
#' @param data A data frame or data table with shape LFDT where all replicated measurements are enclosed
#'
#' @return A data table with format LFDT representing precision estimates for each MS comparison. MS variances, MS CVs and lambda is returned.
#' @export
#'
#' @examples MS_precision(MS_wise(sampled_cs_measurements)[Comparison=="MP1 - MP2",])

MS_precision <- function(data){
  data <- setDT(data)
  names <- names(data)[!names(data)%in%c("ReplicateID","MP_A","MP_B")]
  vars <- data[, .(Var_A = var(MP_A), Var_B = var(MP_B), Mean_A = mean(MP_A), Mean_B = mean(MP_B)),by=names]
  anal <- vars[, .(Var_A = mean(Var_A), Var_B = mean(Var_B), Mean_A = mean(Mean_A), Mean_B = mean(Mean_B))]
  return(anal[, .(Var_A = mean(Var_A), Var_B = mean(Var_B), `CV_A (%)` = 100*sqrt(Var_A)/Mean_A, `CV_B (%)` = 100*sqrt(Var_B)/Mean_B,lambda = mean(Var_A)/mean(Var_B))])
}
