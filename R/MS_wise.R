#' Converting standard wide format MS comparison data (WFDT) to standard long format (LFDT)
#'
#' @param data A data frame or data table with m + n columns where m (usually two) of them are ID columns and n of them is numeric columns consisting of measurement results for each measurement system
#' @param force_IDs If you want to type in the ID column names manually, you may specify them here. This is only recommended if the column ID names are far from the standard names. For example, if the Sample ID column is named 'WDXQXX' and Replicate ID column 'RGDHEC', or is in a language very different from Scandinavian languages and English
#'
#' @details If you specify \code{force_IDs} to your ID column names, the output will keep these names. It is then recommended to run \code{fix_column_names_2(data = output, n_ids = 3, LFDT = TRUE, MOR = FALSE, slack = 4)} to make the output into the standard form so that further functions may be used
#'
#' @return A data table of format LFDT, with m + 3 columns where the latter three columns are MS comparison ID column (new) and resulting numeric columns in this comparison
#' @export
#'
#' @examples MS_wise(sampled_eqam_measurements)

MS_wise <- function(data, force_IDs = NULL){
  IDs <- c("SampleID", "ReplicateID")
  if(!is.null(force_IDs)){
    IDs <- force_IDs
  }
  if(!any(IDs %in% names(data))){
    return(MS_wise_2(data = data, fix_names = TRUE))
  }
  else{
    alphabetically <- c(IDs, sort(setdiff(x = colnames(data), y = IDs)))
    data <- as.data.table(data)[,..alphabetically]
    nord <- c("Comparison", IDs, "MP_A", "MP_B")
    mdt <- melt(data, id.vars = IDs)
    res <- mdt[mdt, on = c(IDs, "variable < variable"), nomatch = NULL,
               c(.(Comparison = paste(x.variable, i.variable, sep = " - "),
                   MP_A = x.value, MP_B = i.value), .SD),
               .SDcols = IDs][,..nord]
    return(res)
  }

}
