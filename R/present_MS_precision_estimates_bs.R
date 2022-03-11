#' Present precision estimates with bootstrapped confidence interval for all unique measurement systems
#'
#' @param data A data frame or data table with format LFDT or WFDT enclosed with all replicated measurements
#' @param groups Character vector - The names of the grouping columns of the data. Default is 'Comparison'
#' @param LFDT Logical - Is data of format LDFT or WFDT (default)?
#' @param decimals Integer - how many decimals should be used in the output data table. Default is 6
#' @param n The number of bootstrap resamples of MS precision values. Note that choosing a large number may result in a significant computation time
#' @param level A numeric value that captures the overall confidence level of bootstrapped confidence intervals
#' @param seed_id Integer - Set a unique seed for reproducibility, or generate a random integer based on a distribution
#'
#' @return A data table with n rows and 3 columns, where n is the number of unique measurement systems in data
#' @export
#'
#' @examples present_MS_precision_estimates_bs(sampled_cs_measurements, decimals = 3, n = 50, level = 0.99, seed_id = round(runif(1,1,100)))
present_MS_precision_estimates_bs <- function(data, groups = "Comparison", LFDT = FALSE, decimals = 6, n = 1e2, level = 0.95, seed_id = 1){
  set.seed(seed_id)
  alpha <- (1 - level)/2
  if(!LFDT){
    data <- MS_wise(data)
  }
  cpr <- resample_MS_precision_over_groups(data = data, n = n, groups = groups)
  setnames(x = cpr, old = c("CV_A (%)", "CV_B (%)"), new = c("CV_A","CV_B"), skip_absent = TRUE)
  cpr <- melt(data = cpr,
             id.vars = groups,
             variable.name = "variable",
             value.name = "value")[variable!="lambda",]
  cpr <- cpr[,.(`Precision measure` = unlist(lapply(X=stri_split(str=variable,fixed="_"),function(x)x[1])),
              MPID = ifelse(stri_sub(variable,-1)=="A",
                            unlist(lapply(X=stri_split(str=Comparison,fixed=" - "),function(x)x[1])),
                            unlist(lapply(X=stri_split(str=Comparison,fixed=" - "),function(x)x[2]))),
              estimate = value)]
  cpr <- cpr[,.(mean = round(mean(estimate), decimals),
                lwr = round(quantile(x = estimate, probs = alpha), decimals),
                upr = round(quantile(x = estimate, probs = 1 - alpha), decimals)),
                by = list(`Precision measure`, MPID)]
  cpr <- cpr[,.(`estimate + confidence interval` = paste0(mean," (",lwr,", ",upr,")")),
             by = list(`Precision measure`, MPID)]
  cpr <- dcast(data = cpr, formula = MPID ~ `Precision measure`, value.var = "estimate + confidence interval")
  setnames(x = cpr, old = c("MPID","CV","Var"), new = c("Name of measurement system", paste0("CV estimate (%) with ",level*100,"% confidence interval"), paste0("Variance estimate with ",level*100,"% confidence interval")))
  cpr
}


