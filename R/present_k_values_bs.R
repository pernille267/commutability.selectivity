#' Present calculated k values with bootstrapped confidence interval for all unique measurement system comparisons
#'
#' @param data A data frame or data table with format LFDT or WFDT (default) enclosed with all replicated measurements
#' @param groups Character vector - The names of the grouping columns of the data. Default is 'Comparison'
#' @param LFDT Logical - Is data of format LDFT or WFDT (default)?
#' @param decimals Integer - how many decimals should be used in the output data table. Default is 6
#' @param n The number of bootstrap resamples of k values for each group. Note that choosing a large number may result in a significant computation time. n > 500 is not recommended!
#' @param level A numeric value that captures the overall confidence level of bootstrapped confidence intervals
#' @param seed_id Integer - Set a unique seed for reproducibility, or generate a random integer based on a distribution
#' @param two_sided Logical - should two-sided confidence intervals be calculated. Default is no because of the practical use of k
#'
#' @return A data table with n rows and 2 columns, where n is the number of unique measurement system combinations in data
#' @export
#'
#' @examples present_k_values_bs(data = sampled_cs_measurements, n = 10, level = 0.99, seed_id = 5, two_sided = TRUE)
present_k_values_bs <- function(data, groups = "Comparison", LFDT = FALSE, decimals = 6, n = 1e2, level = 0.95, seed_id = 1, two_sided = FALSE){
  set.seed(seed_id)
  alpha <- ifelse(test = two_sided, (1 - level) / 2, 1 - level)
  if(!LFDT){
    data <- MS_wise(data)
  }
  ckr <- resample_k_over_groups(data = data, n = n, groups = groups)
  ckr <- ckr[,.(k = round(mean(k),decimals),
         lwr = ifelse(!two_sided,0,round(quantile(x=k,probs=alpha),decimals)),
         upr = round(quantile(x=k,probs=1-alpha),decimals)), by = groups]
  ckr <- ckr[,.(k = paste0(k," (",ifelse(!two_sided,"<-",lwr),", ",upr,")")), by = groups]
  setnames(ckr,"k",paste0("k value with ",level*100,"% ", ifelse(!two_sided,"upper ","two-sided "),"confidence interval"))
  ckr
}

