#' Attach calculated k values for each group
#'
#' @param data A grouped data frame or data table with format LFDT that is the one we want the k values attached to
#' @param k_data A grouped (same grouping variables as data) data frame or data table with format LFDT containing the k values
#' @param groups The names of the grouping columns of data and k_data
#'
#' @return
#' @export
#'
#' @examples attach_k_values_over_groups(data = Deming_CLSI_EP14_over_groups(MS_wise(sampled_cs_measurements), groups = "Comparison", level = 0.95, Np = 1, evaluated_materials = MS_wise(sampled_eqam_measurements)), k_data = estimate_k_over_groups(MS_wise(sampled_cs_measurements)))
attach_k_values_over_groups <- function(data, k_data, groups = "Comparison"){
  data <- as.data.table(data)
  k_data <- as.data.table(k_data)
  return(merge(x = data, y = k_data, by = groups))
}
