#' Prepare data for constructing histograms of residuals
#'
#' @param data A data table or data frame that is an output form \code{residuals_vs_fitted()}
#' @param exclude_fitted Should the fitted values be discarded? Default is \code{TRUE}
#' @param remove_NA Should NA-values be removed from data. Default is \code{TRUE}
#' @param method Which method should be used to obtain the optimal number of bins, that is, the width of the bins. Valid inputs are \code{FDR}, \code{SNRR}, \code{SF}, \code{RR} and \code{DF}. See details for more information regarding the different approaches
#'
#' @return A data table that have new columns relevant for \code{plot_residual_histograms()}
#' @export
#'
#' @examples print(1)

residual_histogram_data <- function(data, exclude_fitted = TRUE, remove_NA = TRUE, method = c("FDR", "SNRR", "SF", "RR", "DF")){
  data <- setDT(data)
  if(exclude_fitted){
    data <- data[,c("Comparison","residuals")]
  }
  data <- data[,list(residuals = residuals,
                     studentized_residuals = studentize(residuals),
                     x = seq(from = min(residuals, na.rm = remove_NA),
                             to = max(residuals, na.rm = remove_NA),
                             length.out = ifelse(remove_NA,
                                                 length(na.omit(residuals)),
                                                 length(residuals))),
                     x_studentized = seq(from = min(studentize(residuals), na.rm = remove_NA),
                                         to = max(studentize(residuals), na.rm = remove_NA),
                                         length.out = ifelse(remove_NA,
                                                             length(na.omit(studentize(residuals))),
                                                             length(studentize(residuals)))),
                     mean_residuals = mean(residuals, na.rm = remove_NA),
                     sd_residuals = sd(residuals, na.rm = remove_NA),
                     optimal_binwidth_raw = get_optimal_binwidth(data = residuals, method = method),
                     optimal_binwidth_studentized = get_optimal_binwidth(data = studentize(residuals), method = method)),
               by = list(Comparison)]
  return(data)
}
