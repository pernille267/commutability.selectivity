#' Automatically fix column names using the second algorithm
#'
#' @param data A data frame or data table with format LFDT or WFDT having \code{n_ids} ID columns. The ID columns are not required to be ordered as is the case for \code{fix_column_names()}
#' @param n_ids The number of ID columns found in \code{data}
#' @param LFDT Is the data on long-format with comparison, i.e., TRUE or on wide format without comparison which corresponds to FALSE
#' @param MOR Is the data only consisting of mean of replicates? If so, set this to TRUE. Or else, set this to FALSE
#' @param slack How much deviance from look-up table should be allow.
#' @param print_matches Should the matches for the respective ID columns from \code{data} be displayed? If there are more than four matches, only the first four will be presented and the rest omitted
#'
#' @description In order to use many of the function of this package, it is required that the ID columns are given by particular names. Comparison ID column should be named \code{Comparison}, Clinical sample ID should be named \code{SampleID} and Replicate measurement ID should be named \code{ReplicateID}. If given names deviates from the standard names, one might experience problems when using some of the functions found in this package
#'
#' @details You need to be careful with how you specify LFDT, MOR and n_ids as wrong combinations of these may not be supported. Contact Author for further information
#'
#' @return A data table identical to input \code{data} but ID column names are potentially replaced with standard ones if they are found in respective look-up tables
#' @export
#'
#' @examples fix_column_names_2(sampled_eqam_measurements, LFDT = FALSE, MOR = FALSE, slack = 1, print_matches = TRUE)

fix_column_names_2 <- function(data, n_ids = 3, LFDT = FALSE, MOR = FALSE, slack = 3, print_matches = FALSE){
  data <- setDT(data)
  y <- names(data)[1:n_ids]
  x <- stri_replace_all(tolower(y),replacement="",fixed = " ")
  if(all(n_ids == 3, LFDT, !MOR)){
    cid <- 1
    sid <- 2
    rid <- 3
  }
  else if(all(n_ids == 2, LFDT, MOR)){
    cid <- 1
    sid <- 2
  }
  else if(all(n_ids == 2, !LFDT, !MOR)){
    sid <- 1
    rid <- 2
  }
  else if(all(n_ids == 1, !LFDT, MOR)){
    sid <- 1
  }
  else{
    stop("no supported combination of n_ids, LFDT and MOR")
  }
  accepted_comparison_id_names <- c("comparison","comp","cmparison","coparison","compirison","comprison","sammenligning","mscomparison","mcomparison","metodesammenligning","mssammenligning")
  accepted_sample_id_names <- c("sampleid","csid","clinicalsampleid","patient","cs","cid","sample","clinid","patientid","sampl","sam","pasient","prøve","prøveid","pasientid")
  accepted_replicate_id_names <- c("replicateid","rid","replicates","replicate","replicatesid","rep","repl","replic","repeatid","repid","replicat","replikatid","replikat","vid")

  sufficiently_close_matches <- list("Comparison" = NA, "SampleID" = NA, "ReplicateID" = NA)
  matching_messages <- list("Comparison match(es) with lookup:" = "No good matches for comparison id column",
                            "SampleID match(es) with lookup:" = "No good matches for sample id column",
                            "ReplicateID match(es) with lookup:" = "No good matches for replicate id column")
  closest_sampleid <- grab(x = accepted_sample_id_names, pattern = x[sid], maxDist = slack, value = TRUE)

  if(length(closest_sampleid)>0){
    if(length(closest_sampleid)<=3){
      matching_messages[[2]] <- stri_c(closest_sampleid, collapse = ", ")
    }
    else if(length(closest_sampleid)>3){
      closest_sampleid <- closest_sampleid[1:4]
      matching_messages[[2]] <- stri_c(c(closest_sampleid, "..."), collapse = ", ")
    }

    sufficiently_close_matches[[2]] <- TRUE
  }
  else{
    sufficiently_close_matches[[2]] <- NA
  }
  if(LFDT){
    closest_comparison <- grab(x = accepted_comparison_id_names, pattern = x[cid], maxDist = slack, value = TRUE)
    if(length(closest_comparison)>0){
      if(length(closest_comparison) <= 3){
        matching_messages[[1]] <- stri_c(closest_comparison, collapse = ", ")
      }
      else{
        closest_comparison <- closest_comparison[1:4]
        matching_messages[[1]] <- stri_c(c(closest_comparison, "..."), collapse = ", ")
      }
      sufficiently_close_matches[[1]] <- TRUE
    }
    else{
      sufficiently_close_matches[[1]] <- NA
    }
  }

  if(!MOR){
    closest_replicateid <- grab(x = accepted_replicate_id_names, pattern = x[rid], maxDist = slack, value = TRUE)
    if(length(closest_replicateid)>0){
      if(length(closest_replicateid) <= 3){
        matching_messages[[3]] <- stri_c(closest_replicateid, collapse = ", ")
      }
      else{
        closest_replicateid <- closest_replicateid[1:4]
        matching_messages[[3]] <- stri_c(c(closest_replicateid,"..."), collapse = ", ")
      }
      sufficiently_close_matches[[3]] <- TRUE
    }
    else{
      sufficiently_close_matches[[3]] <- NA
    }
  }
  if(print_matches){
    matching_messages <- unlist(matching_messages)
    cat(matching_messages, sep = "\n", fill = TRUE, labels = c("Comparison ID column matches the following from LT:  ",
                                                               "Sample ID column matches the following from LT:  ",
                                                               "Replicate ID column matches the following from LT:  "))
  }

  valid_names <- stri_replace_na(str = unlist(sufficiently_close_matches), replacement = FALSE)
  names_to_be_included <- names(sufficiently_close_matches)[which(valid_names==TRUE)]
  names(data)[1:n_ids] <- names_to_be_included
  return(data)
}

