#' Resample data with replicated measurements
#'
#' @param data A data frame or data table with format LFDT where all replicate measurements are enclosed. Two of the columns must be the standard ID columns SampleID and ReplicateID
#'
#' @return A data table that is the resulting resampled data with CSs as clusters
#' @export
#'
#' @examples resample_replicated_data(data = MS_wise(sampled_cs_measurements)[Comparison=="MP1 - MP3",])

resample_replicated_data <- function(data){
  SampleIDs <- unique(data$SampleID)
  ReplicateIDs <- unique(data$ReplicateID)
  n <- length(SampleIDs)
  R <- length(ReplicateIDs)
  rcss <- data.table(Resampled_SampleID = rep(sample(x = SampleIDs, size = n, replace = TRUE), each = R), ReplicateID = rep(ReplicateIDs, times = n))
  data <- merge(data, rcss, by.x = c("SampleID","ReplicateID"), by.y = c("Resampled_SampleID","ReplicateID"))
  data <- data[,`:=`(ReplicateID = rep(ReplicateIDs,times=n))]
  return(data)
}


