#' Calculate k values over multiple groups
#'
#' @param data A grouped data frame or data table with format LFDT
#' @param groups The names of the grouping columns of the data
#'
#' @return A data table listing the calculated k values for each group
#' @export
#'
#' @examples estimate_k_over_groups(data = MS_wise(sampled_cs_measurements), groups = "Comparison")

estimate_k_over_groups <- function(data, groups = "Comparison"){
  data <- as.data.table(data)
  rbindlist(l = lapply(X = split(data, by = groups), FUN = estimate_k), idcol = paste(groups,sep=" - "))
}

