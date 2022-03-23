#' Converting standard wide format MS comparison data (WFDT) to standard long format (LFDT)
#'
#' @param data A data frame or data table with m + n columns where m (usually two) of them are ID columns and n of them is numeric columns consisting of measurement results for each measurement system. The ID columns may be named approximately e.g., SampleID or ReplicateID, but names deviating much from these will not be recongnized
#' @param fix_names Should ID column names with inappropriate names be attempted fixed. Default is \code{TRUE}, if we set \code{FALSE} and the ID columns not named correctly, this function will produce an error
#'
#' @description A backup function or convenience function where we do not need to specify the ID columns explicitly. The given ID column names must however be close to corresponding standard names
#'
#' @return A data table of format LFDT, with m + 3 columns where the latter three columns are MS comparison ID column (new), Sample ID column and Replicate ID column. The m remaining columns are the resulting numeric columns representing MS measurements
#' @export
#'
#' @examples MS_wise_2(sampled_cs_measurements, fix_names = FALSE)

MS_wise_2 <- function(data, fix_names = TRUE){
  data <- setDT(data)
  IDs <- c("SampleID", "ReplicateID")
  if(fix_names){
    data <- fix_column_names_2(data = data,
                               n_ids= 2,
                               LFDT = FALSE,
                               MOR  = FALSE,
                               slack= 4,
                               print_matches = FALSE)
    IDs <- names(data)[which(names(data)%in%IDs)]
    if(length(IDs)==0){
      stop("No close approximate matches for ID columns")
    }
  }
  alphabetically <- c(IDs, sort(setdiff(x = names(data), y = IDs)))
  data <- data[,..alphabetically]
  nord <- c("Comparison", IDs, "MP_A", "MP_B")
  mdt <- melt(data, id.vars = IDs)
  res <- mdt[mdt, on = c(IDs, "variable < variable"), nomatch = NULL,
             c(.(Comparison = paste(x.variable, i.variable, sep = " - "),
                 MP_A = x.value, MP_B = i.value), .SD),
             .SDcols = IDs][,..nord]
  return(res)
}
