#' resample values in one specific group
#'
#' @param data a subset of the total data
#' @param B number of bootstrap replicates
#' @param fun which statistic
#'
#' @return return something cool
#' @export
#'
#' @examples print(1)
gro_resample <- function(data, B = 1000, fun = median){
  new_data <- pbreplicate(n = B, expr = data[,.(RD=sample(RD,length(RD),TRUE))],simplify=FALSE)
  new_data <- rbindlist(l = new_data, idcol = "RSTEP")
  return(new_data)
}

