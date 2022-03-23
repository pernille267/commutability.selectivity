#' Automatically fix column names using the first algorithm
#'
#' @param data A data frame or data table with format LFDT or WFDT having \code{n_ids} ID columns
#' @param LFDT Is the data on long-format with comparison, i.e., TRUE or on wide format without comparison which corresponds to FALSE
#' @param MOR Is the data only consisting of mean of replicates? If so, set this to TRUE. Or else, set this to FALSE
#' @param ordered Is the data ordered with \code{Comparison} first, then \code{SampleID}, and last \code{ReplicateID}? If not, data is sent to \code{fix_column_names_2()}, which does approximate matching with slack = 3. Note that setting to ordered when not being ordered in the correct way will affect the output data in some negative manner
#' @param n_ids The number of ID columns found in \code{data}. Ignored if \code{ordered} is \code{TRUE}
#'
#' @return A data table identical to input \code{data} but ID column names are replaced with standard ones. If ordered is \code{FALSE} then approximate matching is done with slack = 3
#' @export
#'
#' @examples fix_column_names(data = sampled_cs_measurements, ordered = TRUE)

fix_column_names <- function(data, LFDT = FALSE, MOR = FALSE, ordered = TRUE, n_ids = 3){
  if(ordered){
    if(all(LFDT,!MOR)){
      names(data)[1:3] <- c("Comparison","SampleID","ReplicateID")
    }
    else if(all(!LFDT,!MOR)){
      names(data)[1:2] <- c("SampleID", "ReplicateID")
    }
    else if(all(LFDT,MOR)){
      names(data)[1:2] <- c("Comparison", "SampleID")
    }
    else if(all(!LFDT,MOR)){
      names(data)[1] <- c("SampleID")
    }
    else{
      stop("no supported combination of LFDT and MOR")
    }
    return(data)
  }
  else{
    return(fix_column_names_2(data = data, n_ids = n_ids, LFDT = LFDT, MOR = MOR, slack = 3, print_matches = FALSE))
  }

}

