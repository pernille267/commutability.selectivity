#' Calculate mean of replicates over Sample IDs and possibly other grouping columns
#'
#' @param data A data frame or data table with format LFDT where all replicated measurements are enclosed. If \code{Data} is on WFDT, you must use \code{mean_of_replicates_2()} instead
#'
#' @description Calculate mean of replicates based on \code{data}. Note that the three columns \code{ReplicateID}, \code{MP_A} and \code{MP_B} must be part of \code{data} for this function to work properly. The grouping columns may have any name, but it is recommended to stick to standard names such as SampleID or Comparison
#'
#' @return A data table with format LFDT where replicated measurements are collapsed to mean of replicated measurements
#' @export
#'
#' @examples mean_of_replicates(data = MS_wise(sampled_eqam_measurements)[Comparison=="MP1 - MP2",])

mean_of_replicates <- function(data){
  data <- setDT(data)
  names <- names(data)[!names(data) %in% c("ReplicateID", "MP_A", "MP_B")]
  data <- data[, list(MP_A = mean(MP_A), MP_B = mean(MP_B)), by = names]
  return(data)
}




