#' Calculate mean of replicates based on grouping
#'
#' @param data A data frame or data table with format LFDT where all replicated measurements are enclosed
#' @param IDs The names of the grouping columns of the data. The column of replicate IDs must not be part of these names
#'
#' @return A data table with format LFDT where replicated measurements are collapsed to mean of replicated measurements
#' @export
#'
#' @examples mean_of_replicates(data = MS_wise(sampled_eqam_measurements)[Comparison=="MP1 - MP2",])

mean_of_replicates <- function(data,IDs=c("Comparison","SampleID")){
  data <- as.data.table(data)
  names <- names(data)[names(data)!="ReplicateID"]
  sdcols <- names[!names %in% IDs]
  data <- data[,..names][,lapply(X=.SD,FUN=mean),.SDcols=sdcols,by=IDs]
  return(data)
}


