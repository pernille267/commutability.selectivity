#' Calculate mean of replicates over Sample IDs and possibly other grouping columns specified
#'
#' @param data A data frame or data table with format WFDT where all replicated measurements are enclosed
#' @param other_ids If additional grouping columns should be included, you include them by name here. If your SampleID and ReplicateID columns are named something else, then you need to specify their actual names here
#'
#' @return A data table with format WFDT where replicated measurements are collapsed to mean of replicated measurements
#' @export
#'
#' @examples mean_of_replicates_2(sampled_cs_measurements)

mean_of_replicates_2 <- function(data, other_ids = NULL){
  data <- setDT(data)
  if(is.null(other_ids)){
    MS_columns <- names(data)[!names(data) %in% c("SampleID", "ReplicateID")]
    data <- data[, lapply(X = .SD, fun = mean), .SDcols = MS_columns, by = c("SampleID")]
    return(data)
  }
  else{
    MS_columns <- names(data)[!names(data) %in% c("SampleID", "ReplicateID",other_ids)]
    data <- data[, lapply(X = .SD, fun = mean), .SDcols = MS_columns, by = c("SampleID", other_ids)]
    return(data)
  }
}
