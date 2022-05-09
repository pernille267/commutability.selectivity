#' Estimate prediction band using Deming regression formulated by J. Gillard and Fuller over groups
#'
#' @param data A data frame or data table with format LFDT enclosed with all replicated measurements
#' @param groups The names of the grouping columns of the data
#' @param level A numeric value that captures the overall confidence level of the estimated prediction band
#' @param R An integer signifying the maximum number of replicated measurements performed on each evaluated material
#' @param Np An integer, which captures the number of pointwise prediction intervals making the prediction band across the concentration range. Not relevant if evaluated_materials are !NULL
#' @param evaluated_materials A data frame or data table with format LFDT enclosed with all replicated measurements for evaluated materials such as EQAMs or CRMs. Should be NULL if the PB, that is, pointwise prediction intervals to be estimated
#' @param column_order A vector specifying the order of the columns of the outputted grouped data table. Will be ignored if evaluated_materials is NULL. Note that the length of column_order must be the same as the number of columns of output columns, which is 6 + the number of grouping columns
#' @param add_judgement Should we include a column to the output column stating whether each evaluated material at in each group is inside the estimated prediction interval? Default is \code{FALSE}
#' @param jugdement_labels A character vector including the labels for when EQAM is outside estimated PI (\code{FALSE}) or inside (\code{TRUE}) in that order. For example, judgement_labels = c("Not inside","inside") or judgement_labels = c("No","Yes"). Default is c("no","yes)
#'
#' @return A grouped data table enclosed with information regarding the estimated prediction band across the concentration range or for the particular evaluated materials for all groups
#' @export
#'
#' @details This is the prediction interval estimation procedure using parts from Jonathan Gillard's PhD thesis work and Wayne Fuller's book Measurement error models. This is recommended to use over Deming_CLSI_EP14_over_groups if k is most probably larger than 1
#'
#' @examples Deming_Gillard_Fuller_over_groups(MS_wise(sampled_cs_measurements), groups = "Comparison", level = 0.95, Np = 1, evaluated_materials = MS_wise(sampled_eqam_measurements))

Deming_Gillard_Fuller_over_groups <- function(data, groups = "Comparison", level = 0.99, R = 3, Np = 1e3, evaluated_materials = NULL, column_order = c("Comparison", "SampleID", "MP_B", "MP_A", "fit", "lwr" , "upr"), add_judgement = FALSE, judgement_labels = c("no","yes")){
  if(!is.data.table(data)){
    data <- as.data.table(data)
  }
  data_list <- split(x = data, by = groups)
  if(!is.null(evaluated_materials)){
    if(!any(class(evaluated_materials)[1] == c("data.frame","data.table","tibble","matrix","array","list"))){
      message("For maintainer (evaluated_materials): evaluated_materials is not of correct class")
      stop("evaluated_materials is not data frame or data table. Please make sure that data's class is one of the two")
    }
    if(class(evaluated_materials)[1] == "matrix" | class(evaluated_materials)[1] == "array" | class(evaluated_materials)[1] == "list"){
      if(is.null(colnames(data))){
        message("For maintainer (data): column names of a matrix / array / list must not be NULL !")
        stop("Column names of a matrix / array / list must not be NULL! Make sure they are named by given standards of LFDT")
      }
    }
    evaluated_materials <- as.data.table(evaluated_materials)
    names_evaluated_materials <- names(evaluated_materials)
    names_clinical_samples <- names(data)

    evaluated_materials <- split(x = evaluated_materials, by = groups, keep.by = FALSE)
    if(!all(names_evaluated_materials %in% names_clinical_samples)){
      message("For maintainer: Some names in evaluated_materials does not exist in data")
      stop("evaluated_materials must have identical MS comparisons as data")
    }
    else if(!all(names_evaluated_materials == names_clinical_samples)){
      warning("The order of the columns in evaluated_materials and data is not identical, results may not be trustworthy")
    }
    pb_at_eqmas <- mapply(FUN = function(clinical_samples, evaluated_samples) Deming_Gillard_Fuller(data = clinical_samples, evaluated_materials = evaluated_samples),
                          data_list,
                          evaluated_materials,
                          SIMPLIFY = FALSE)
    pb_at_eqmas <- rbindlist(l = pb_at_eqmas,
                             idcol = paste(groups, sep = " - "))
    setnames(x = pb_at_eqmas,
             old = "X.GRID",
             new = "MP_B",
             skip_absent = TRUE)
    evaluated_materials <- rbindlist(l = lapply(X = evaluated_materials,
                                                FUN = mean_of_replicates),
                                     idcol = paste(groups, sep = " - "))

    output <- merge(x = pb_at_eqmas,
                    y = evaluated_materials,
                    by = c(groups, "MP_B"))

    setorder(x = output, Comparison, SampleID)
    if(ncol(output) == length(column_order)){
      if(all(names(output) %in% column_order)){
        output <- output[,..column_order]
      }
      else{
        message("For maintainer: The names found in output does not all match with the elements of column_order")
        stop("The names found in output does not all match with the elements of column_order. Call is aborted")
      }
    }
    if(add_judgement){
      `inside estimated PI` <- ifelse(output$MP_A > output$lwr & output$MP_A < output$upr, judgement_labels[2], judgement_labels[1])
      output <- cbind(output, "inside estimated PI" = `inside estimated PI`)
    }
    return(unique(output))
  }
  else{
    rbindlist(l = lapply(X = data_list, FUN = Deming_Gillard_Fuller, level = level, R = R, Np = Np, evaluated_materials = NULL), idcol = paste(groups, sep = " - "))
  }
}

