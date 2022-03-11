#' Present precision estimates for unique measurement systems
#'
#' @param data A data frame or data table with format LFDT or WFDT enclosed with all replicated measurements
#' @param groups Character vector - The names of the grouping columns of the data. Default is 'Comparison'
#' @param LFDT Logical - Is data of format LDFT or WFDT (default)?
#' @param decimals Integer - how many decimals should be used in the output data table. Default is 6
#'
#' @return A data table with two rows and 1 + n columns, where n is the number of unique measurement systems in data
#' @export
#'
#' @examples present_MS_precision_estimates(sampled_cs_measurements)
present_MS_precision_estimates <- function(data, groups = "Comparison", LFDT = FALSE, decimals = 6){
  if(!LFDT){
    data <- MS_wise(data = data)
  }
  cp <- MS_precision_over_groups(data = data, groups = groups)
  setnames(x = cp, old = c("CV_A (%)", "CV_B (%)"), new = c("CV_A","CV_B"), skip_absent = TRUE)
  cp <- melt(data = cp,
             id.vars = groups,
             variable.name = "variable",
             value.name = "value")[variable!="lambda",]
  cp <- cp[,.(`Precision measure ` = unlist(lapply(X=stri_split(str=variable,fixed="_"),function(x)x[1])),
              MPID = ifelse(stri_sub(variable,-1)=="A",
                           unlist(lapply(X=stri_split(str=Comparison,fixed=" - "),function(x)x[1])),
                           unlist(lapply(X=stri_split(str=Comparison,fixed=" - "),function(x)x[2]))),
              estimate = value)]
  cp <- unique(cp)[,lapply(X=.SD,FUN=round,decimals=decimals),.SDcols="estimate",by=list(`Precision measure`, `MPID`)]
  cp <- dcast(data = cp, formula = `Precision measure ` ~ MPID, value.var = "estimate")
  cp$`Precision measure ` <- c("MS CV(s) in %:", "MS variance(s):")
  return(cp)
}




