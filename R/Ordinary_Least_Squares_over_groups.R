#' Estimate prediction band or prediction intervals using ordinary least squares regression over groups
#'
#' @param data A grouped data frame or data table with format LFDT enclosed with all replicated measurements
#' @param groups The names of the grouping columns of \code{data} and \code{evaluated_materials}
#' @param level A numeric value that captures the overall confidence level of the estimated prediction band
#' @param R An integer signifying the maximum number of replicated measurements performed on each evaluated material
#' @param Np An integer, which captures the number of pointwise prediction intervals making the prediction band across the concentration range. Not relevant if evaluated_materials are !NULL
#' @param evaluated_materials A data frame or data table with format LFDT enclosed with all replicated measurements for evaluated materials such as EQAMs or CRMs. Should be NULL if the PB, that is, pointwise prediction intervals to be estimated
#' @param optimize_for_measurement_errors Logical - Should we shift axes if the measurement error for MP_B is larger than MP_A? This will be optimal because OLS assumes that variability in MP_B is zero
#' @param column_order A vector specifying the order of the columns of the outputted grouped data table. Will be ignored if evaluated_materials is NULL. Note that the length of column_order must be the same as the number of columns of output columns, which is 6 + the number of grouping columns
#'
#' @details This function is not finished
#'
#' @return A grouped data table enclosed with information regarding the estimated prediction band across the concentration range or prediction intervals for the particular evaluated materials for all groups
#' @export
#'
#' @examples Ordinary_Least_Squares_over_groups(MS_wise(sampled_cs_measurements), level = 0.99, evaluated_materials = MS_wise(sampled_eqam_measurements))
Ordinary_Least_Squares_over_groups <- function(data, groups = "Comparison", level = 0.99, R = 3, Np = 1e3, evaluated_materials, optimize_for_measurement_errors = TRUE, column_order = c("reversed","Comparison", "SampleID", "MP_B", "MP_A", "fit", "lwr" , "upr")){
  data <- split(x = as.data.table(data), by = groups)
  if(!is.null(evaluated_materials)){
    evaluated_materials <- split(x = as.data.table(evaluated_materials), by = groups)
    if(!all(names(evaluated_materials)%in%names(data))){stop("evaluated_materials must have identical MS comparisons as data")}
  }
  if(!is.null(evaluated_materials)){
    pb_at_eqmas <- mapply(FUN = function(x,y) Ordinary_Least_Squares(data = x, evaluated_materials = y, optimize_for_measurement_errors = optimize_for_measurement_errors),
                          data,
                          evaluated_materials,
                          SIMPLIFY = FALSE)
    pb_at_eqmas <- rbindlist(l = pb_at_eqmas, idcol = paste(groups, sep = " - "))
    setnames(x = pb_at_eqmas, old = "X.GRID", new = "MP_B", skip_absent = TRUE)
    evaluated_materials <- rbindlist(l = lapply(X = evaluated_materials,
                                                FUN = function(x) mean_of_replicates(x)[,-1]),
                                                idcol = paste(groups, sep = " - "))
    if(all(optimize_for_measurement_errors, any(pb_at_eqmas$reversed))){
      pb_at_eqams_list <- split(x = pb_at_eqmas, by = "reversed")
      setnames(x = pb_at_eqams_list$`TRUE`, old = "MP_B", new = "MP_A", skip_absent = TRUE)
      output <- merge(x = pb_at_eqmas, y = evaluated_materials, by = c(groups, "MP_B"))
      reversed_part <- merge(x = pb_at_eqams_list$`TRUE`, y = evaluated_materials, by = c(groups, "MP_A"))
      column_order_default <- names(output)
      reversed_part <- reversed_part[,..column_order_default]
      if(length(pb_at_eqams_list)>=2){
        other_part <- merge(x = pb_at_eqams_list$`FALSE`, y = evaluated_materials, by = c(groups, "MP_B"))
        output <- rbind(reversed_part, other_part)
      }
      else{
        output <- reversed_part
      }

    }
    else{
      output <- merge(x = pb_at_eqmas, y = evaluated_materials, by = c(groups, "MP_B"))
    }
    setorderv(x = output, c(groups, "SampleID"), order = rep(1, length(c(groups, "SampleID"))))
    if(ncol(output) == length(column_order)){output <- output[,..column_order]}
    return(unique(fix_reverse_issues(output)))
  }
  else{
    output <- rbindlist(l = lapply(X = data, FUN = Ordinary_Least_Squares, level = level, R = R, Np = Np, evaluated_materials = NULL, optimize_for_measurement_errors = optimize_for_measurement_errors), idcol = paste(groups, sep = " - "))
    output
  }
}
