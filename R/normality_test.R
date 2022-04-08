#' Checks normality of error terms by using Shapiro Wilk's test on each groups' residuals
#'
#' @param data A data table or data frame of type LFDT grouped according to \code{groups}
#' @param alpha The base significance level of the test.
#' @param groups The names of grouping columns in \code{data}
#' @param bonferroni_correct Should we Bonferroni correct the significance level based the number of MS comparisons. Default is \code{TRUE}
#' @param conservative Should the Bonferroni correction be done conservatively that assumes that some of the MS comparisons are dependent? Default is \code{TRUE}
#'
#' @return A data table with grouping as the first column, LM test statistic as the second column, p-value as the third column, and judgement as the last column
#' @export
#'
#' @examples print(1)

normality_test <- function(data, alpha = 0.05, groups = "Comparison", bonferroni_correct = FALSE, conservative = TRUE){
  data <- setDT(data)
  data <- residual_vs_fitted(data = data, groups = groups, MOR = TRUE)
  data_list <- split(x = data, by = groups)
  column_names <- c("statistic", "p.value")
  m <- length(data_list)
  if(bonferroni_correct){
    if(conservative){
      alpha <- alpha / sqrt(m)
    }
    else{
      alpha <- alpha / m
    }
  }
  test_list <- lapply(X = data_list, FUN = function(x) shapiro.test(x = x$residuals))
  test_list <- lapply(X = test_list, FUN = setDT)
  test_list <- lapply(X = test_list, FUN = function(x) x[,..column_names])
  test_list <- lapply(X = test_list,
                      FUN = function(x) x[,list(statistic = statistic,
                                                p.value = p.value,
                                                `Significance level` = alpha,
                                                Judgement = ifelse(p.value < alpha,
                                                                   "Significant non-normality",
                                                                   "normality"))])
  test_results <- rbindlist(l = test_list, use.names = TRUE, idcol = paste0(groups, collapse = " - "), fill = TRUE)
  setnames(test_results, old = column_names, new = c("Shapiro Wilk Statistic", "p-value"))
  return(test_results)
}

