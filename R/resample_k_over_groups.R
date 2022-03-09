#' Resample k values based on grouped data
#'
#' @param data A grouped data frame or data table with format LFDT
#' @param n The number of bootstrap resamples of k for each group. Note that choosing a large number may result in a significant computation time
#' @param groups The names of the grouping columns of the data
#'
#' @return A grouped data table listing the n resampled k values for each group
#' @export
#'
#' @examples resample_k_over_groups(data = MS_wise(sampled_cs_measurements), n = 10)

resample_k_over_groups <- function(data, n=1e2, groups = "Comparison"){
  data <- as.data.table(data)
  rbindlist(l = lapply(X = split(data, by = groups), FUN = resample_k, n = n), idcol = paste(groups,sep=" - "))
}

