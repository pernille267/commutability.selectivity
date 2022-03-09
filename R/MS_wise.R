#' Converting standard wide format MS comparison data (WFDT) to standard long format (LFDT)
#'
#' @param data A data frame or data table with m + n columns where m (usually two) of them are ID columns and n of them is numeric columns consisting of measurement results for each measurement system
#' @param IDs A character vector of length m containing the ID column names
#'
#' @return A data table of format LFDT, with m + 3 columns where the latter three columns are MS comparison ID column (new) and resulting numeric columns in this comparison
#' @export
#'
#' @examples MS_wise(sampled_eqam_measurements, IDs = c("SampleID","ReplicateID"))

MS_wise <- function(data, IDs = c("SampleID","ReplicateID")){
  alphabetically <- c(IDs, sort(setdiff(x = colnames(data), y = IDs)))
  data <- as.data.table(data)[,..alphabetically]
  nord <- c("Comparison", IDs, "MP_A", "MP_B")
  mdt <- melt(as.data.table(data), id.vars = IDs)
  res <- mdt[mdt, on = c(IDs, "variable < variable"), nomatch = NULL,
             c(.(Comparison = paste(x.variable, i.variable, sep = " - "),
                 MP_A = x.value, MP_B = i.value), .SD),
             .SDcols = IDs][,..nord]
  return(res)
}
