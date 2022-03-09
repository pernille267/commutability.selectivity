#' Resample k based on data
#'
#' @param data A data frame or data table with format LFDT with all replicate measurements enclosed
#' @param n The number of bootstrap resamples of k. Note that choosing a large number may result in a significant computation time
#'
#' @return A data table with all resampled k values based on input data
#' @export
#'
#' @examples resample_k(data = MS_wise(sampled_cs_measurements)[Comparison=="MP4 - MP5",])

resample_k <- function(data,n=1e2){
  data <- as.data.table(data)
  all_replicated_dts <- replicate(n = n, expr = resample_replicated_data(data), simplify = FALSE)
  rbindlist(l=lapply(X = all_replicated_dts, FUN = estimate_k))
}



