#' Summarize to residuals vs fitted for each group
#'
#' @param data A data table or data frame with LFDT format
#' @param groups Column name(s) of grouping columns. Default is Comparison
#' @param MOR Should mean of replicates be applied to data before residuals and fitted values are calculated? Default is \code{TRUE}
#'
#' @return A grouped data table according to \code{groups}, where in each group, residuals and fitted values are listed
#' @export
#'
#' @examples residual_vs_fitted(MS_wise(sampled_cs_measurements))

residual_vs_fitted <- function(data, groups = "Comparison", MOR = TRUE){
  data <- as.data.table(data)
  if(MOR){
    data <- mean_of_replicates(data)
  }
  data_list <- split(x = data, by = groups)
  model_list <- lapply(X = data_list, FUN = function(x) lm.fit(x = cbind(1, x$MP_B), y = x$MP_A))
  model_residuals <- lapply(X = model_list, FUN = function(model_data) model_data$residuals)
  model_fitted <- lapply(X = model_list, FUN = function(model_data) model_data$fitted.values)
  residuals_and_fitted <- mapply(function(x,y) setDT(list("residuals" = x, "fitted values" = y)), model_residuals, model_fitted, SIMPLIFY = FALSE)
  residuals_and_fitted <- rbindlist(l = residuals_and_fitted, use.names = TRUE, fill = TRUE, idcol = paste0(groups, collapse = " - "))
  return(residuals_and_fitted)
}

