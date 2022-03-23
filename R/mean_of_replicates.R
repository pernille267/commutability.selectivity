#' Calculate mean of replicates over Sample IDs and possibly other grouping columns if they exist
#'
#' @param data A data frame or data table with format LFDT where all replicated measurements are enclosed. If \code{Data} is on WFDT, you must use \code{mean_of_replicates_2()} instead
#'
#' @description Calculate mean of replicates based on \code{data}. Note that the three columns \code{ReplicateID}, \code{MP_A} and \code{MP_B} must be part of \code{data} for this function to work properly. The grouping columns may have any name, but it is recommended to stick to standard names such as SampleID or Comparison
#'
#' @details If the column displaying replicate IDs are named something else than \code{ReplicateID}, it will be approximately matched, but there is no guarantee that this will work
#'
#' @return A data table with format LFDT where replicated measurements are collapsed to mean of replicated measurements
#' @export
#'
#' @examples mean_of_replicates(data = MS_wise(sampled_eqam_measurements)[Comparison=="MP1 - MP2",])

mean_of_replicates <- function(data){
  data <- setDT(data)
  if(!"ReplicateID"%in%names(data)){
    proposed_id <- which(stringdist::ain(x = names(data), table = c("replicatid","replic","replica","repl","rep","re","repeat","repet","replika","replikat","replikk","replukat","replakat","replicite","replicatt","rreplicateid","repldid","eplicateid","eplicat","feplicatei"), maxDist = 4, method = "lcs"))[1]
    if(length(proposed_id)<1){
      stop("Nothing close to a column with name ReplicateID are found. Please make sure you have one column named like this")
    }
    else{
      names(data)[proposed_id] <- "ReplicateID"
    }
  }
  names <- names(data)[!names(data) %in% c("ReplicateID", "MP_A", "MP_B")]
  data <- data[, list(MP_A = mean(MP_A), MP_B = mean(MP_B)), by = names]
  return(data)
}




