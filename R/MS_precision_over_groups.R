#' Estimate MS precision over multiple groups
#'
#' @param data A grouped data frame or data table with format LFDT
#' @param groups The names of the grouping columns of the data
#'
#' @return A data table listing the estimated MS precision for each group
#' @export
#'
#' @examples MS_precision_over_groups(data = MS_wise(sampled_cs_measurements))

MS_precision_over_groups <- function(data, groups = "Comparison"){
  data <- as.data.table(data)
  rbindlist(l = lapply(X = split(data, by = groups), FUN = MS_precision), idcol = paste(groups,sep=" - "))
}

