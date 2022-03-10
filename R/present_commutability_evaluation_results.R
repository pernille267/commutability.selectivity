#' Present every evaluated materials' commutability for every unique measurement system
#'
#' @param data A data frame or data table with format LFDT or WFDT enclosed with all replicated measurements. This is needed for computation
#' @param evaluated_materials A data frame or data table with format LFDT or WFDT enclosed with all replicated measurements for evaluated materials such as EQAMs or CRMs. This is needed for computation
#' @param LFDT Is data and evaluated_materials of format LDFT or WFDT?
#' @param R An integer signifying the maximum number of replicated measurements performed on each evaluated material. This will shrink the width of the prediction interval and band because using the mean of replicates decreases uncertainty
#' @param level A numeric value that captures the overall confidence level of the estimated prediction intervals. This is assumed to already be Bonferroni-corrected. The recommended base level is 0.99
#' @param method Which regression estimator should be used to calculate prediction intervals. Write "DFG" (Default) for Deming regression with Gillard and Fuller method. Write "DCE" for CLSI method
#' @param success In the summary table, what should we denote those pairs of evaluated materials and measurement systems where the evaluated materials are within prediction interval limits
#' @param failure In the summary table, what should we denote those pairs of evaluated materials and measurement systems where the evaluated materials are outside prediction interval limits
#'
#' @return A data table with 1 + n columns and m rows. Here n is the number of unique measurement systems considered and m is the number of unique evaluated materials
#' @export
#'
#'
#' @examples present_commutability_evaluation_results(data = sampled_cs_measurements, evaluated_materials = sampled_eqam_measurements, method = "DCE")
present_commutability_evaluation_results <- function(data, evaluated_materials, LFDT = FALSE, R = 3, level = 0.99, method = "DGF", success = "Commutable", failure = "Non-commutable"){
  data <- as.data.table(data)
  evaluated_materials <- as.data.table(evaluated_materials)
  csl <- data
  eql <- evaluated_materials
  if(!LFDT){
    csl <- MS_wise(data = data)
    eql <- MS_wise(data = evaluated_materials)
  }
  if(method == "DFG"){
    cel <- Deming_Gillard_Fuller_over_groups(data = csl,
                                             level = level,
                                             R = R, Np = 1,
                                             evaluated_materials = evaluated_materials)
  }
  else if(method != "DFG"){
    cel <- Deming_CLSI_EP14_over_groups(data = csl,
                                        level = level,
                                        R = R, Np = 1,
                                        evaluated_materials = evaluated_materials)
  }
  cel <- cel[,.(inside = MP_A > lwr & MP_A < upr),
             by = list(Comparison, SampleID)]
  cel <- cel[,.(SampleID = SampleID,
                MPID = unlist(stri_split(str = Comparison, fixed = " - ")),
                inside = inside)]
  cel <- cel[,.(`Commutable for particular MS` = ifelse(all(inside), success, failure)),
             by = list(MPID, SampleID)]
  setorderv(x = cel,
            cols = c("SampleID", "MPID"),
            order = c(1,1))
  dcast.data.table(data = cel, formula = SampleID ~ MPID, value.var = "Commutable for particular MS")
}




