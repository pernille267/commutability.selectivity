#' Testing for variance-inhomogeneity using Breauch-Pagan test for each comparison
#'
#' @param data A data table or data frame of type LFDT grouped according to \code{groups}
#' @param alpha The base significance level of the test.
#' @param groups The names of grouping columns in \code{data}
#' @param bonferroni_correct Should we Bonferroni correct the significance level based the number of MS comparisons. Default is \code{TRUE}
#' @param conservative Should the Bonferroni correction be done conservatively that assumes that some of the MS comparisons are dependent? Default is \code{TRUE}
#' @param robust Should Koenker's adjust be implemented in the test? Recommended if data is non-normal. Default is \code{TRUE}
#'
#' @return A data table with grouping as the first column, LM test statistic as the second column, p-value as the third column, and judgement as the last column
#' @export
#'
#' @examples print(1)
heteroscedasticity_test <- function(data, alpha = 0.05, groups = "Comparison", bonferroni_correct = FALSE, conservative = TRUE, robust = TRUE){
  data_list <- split(setDT(data), by = groups)
  data_list <- lapply(X = data_list, FUN = mean_of_replicates)
  m <- length(data_list)
  column_names <- c("statistic", "p.value")
  if(bonferroni_correct){
    if(conservative){
      alpha <- alpha / sqrt(m)
    }
    else{
      alpha <- alpha / m
    }
  }
  test_list <- lapply(X = data_list, FUN = function(x) breusch_pagan(mainlm = list("X" = cbind(1, x$MP_B), "y" = x$MP_A), koenker = robust))
  test_list <- lapply(X = test_list, FUN = setDT)
  test_list <- lapply(X = test_list, FUN = function(x) x[,..column_names])
  test_list <- lapply(X = test_list, FUN = function(x) x[,list(statistic = statistic, p.value = p.value, `Significance level` = alpha, Judgement = ifelse(p.value < alpha, "Significant heteroscedasticity", "homoscedasticity"))])
  test_results <- rbindlist(l = test_list, use.names = TRUE, idcol = paste0(groups, collapse = " - "), fill = TRUE)
  setnames(test_results, old = column_names, new = c("Breusch Pagan Statistic", "p-value"))
  return(test_results)
}

