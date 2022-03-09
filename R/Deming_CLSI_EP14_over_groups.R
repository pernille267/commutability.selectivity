#' Estimate prediction band using Deming regression formulated as in EP14 over groups
#'
#' @param data A data frame or data table with format LFDT enclosed with all replicated measurements
#' @param groups The names of the grouping columns of the data
#' @param level A numeric value that captures the overall confidence level of the estimated prediction band
#' @param R An integer signifying the maximum number of replicated measurements performed on each evaluated material
#' @param Np An integer, which captures the number of pointwise prediction intervals making the prediction band across the concentration range
#' @param evaluated_materials A data frame or data table with format LFDT enclosed with all replicated measurements for evaluated materials such as EQAMs or CRMs. Should be NULL if the PB, that is, pointwise prediction intervals to be estimated
#' @param column_order A vector specifying the order of the columns of the outputted grouped data table. Will be ignored if evaluated_materials is NULL. Note that the length of column_order must be the same as the number of columns of output columns, which is 6 + the number of grouping columns
#'
#' @return A grouped data table enclosed with information regarding the estimated prediction band across the concentration range or for the particular evaluated materials for all groups
#' @export
#'
#' @details This is the prediction interval estimation procedure as described by CLSI in EP14. This is recommended to use over Deming_Gillard_Fuller_over_groups if the calculated k value is smaller than 1
#'
#' @examples Deming_CLSI_EP14_over_groups(MS_wise(sampled_cs_measurements), groups = "Comparison", level = 0.95, Np = 1, evaluated_materials = MS_wise(sampled_eqam_measurements))

Deming_CLSI_EP14_over_groups <- function(data, groups = "Comparison", level = 0.99, R = 3, Np = 1e3, evaluated_materials = NULL, column_order = c("Comparison", "SampleID", "MP_B", "MP_A", "fit", "lwr" , "upr")){
  data <- split(x = as.data.table(data), by = groups)
  if(!is.null(evaluated_materials)){
    evaluated_materials <- split(x = as.data.table(evaluated_materials), by = groups)
    if(!all(names(evaluated_materials) %in% names(data))){stop("evaluated_materials must have identical MS comparisons as data")}
  }
  if(!is.null(evaluated_materials)){
    pb_at_eqmas <- mapply(FUN = function(x,y) Deming_CLSI_EP14(data = x, evaluated_materials = y), data, evaluated_materials, SIMPLIFY = FALSE)
    pb_at_eqmas <- rbindlist(l = pb_at_eqmas, idcol = paste(groups, sep = " - "))
    setnames(x = pb_at_eqmas, old = "X.GRID", new = "MP_B", skip_absent = TRUE)
    evaluated_materials <- rbindlist(l = lapply(X = evaluated_materials, FUN = function(x) mean_of_replicates(x)[,-1]), idcol = paste(groups, sep = " - "))
    output <- merge(x = pb_at_eqmas, y = evaluated_materials, by = c(groups, "MP_B"))
    setorder(x = output, Comparison, SampleID)
    if(ncol(output) == length(column_order)){output <- output[,..column_order]}
    return(unique(output))
  }
  else{
    rbindlist(l = lapply(X = data, FUN = Deming_CLSI_EP14, level = level, R = R, Np = Np, evaluated_materials = NULL), idcol = paste(groups, sep = " - "))
  }

}



