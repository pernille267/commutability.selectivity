#' Estimate precision estimates for each clinical sample over groups
#'
#' @param data A data table or data frame with LFDT format. Must at least contain MP_A, MP_B and something that resembles ID columns
#' @param measure Which precision measure should be used to calculate precision of clinical samples against concentration range? Valid inputs are 'standard deviation' (alternatively, 1 or SD), 'coefficient of variation' (alternatively, 2 or CV) or 'lambda' (alternatively, 3 or L)
#'
#' @return Returns a grouped data table with concentration and precision estimates. If \code{measure = 'lambda'}, the MS column will not be included
#' @export
#'
#' @examples precision_of_replicates(MS_wise(sampled_cs_measurements), measure = 'cv')
precision_of_replicates <- function(data, measure = c("standard deviation","coefficient of variation","lambda")){
  data <- setDT(data)
  data_names <- names(data)
  if(!all(c("SampleID","ReplicateID") %in% data_names)){
    data <- fix_column_names_2(data = data, n_ids = 3, LFDT = TRUE, MOR = FALSE, slack = 4, print_matches = FALSE)
  }
  data_names <- names(data)
  if(!"ReplicateID" %in% data_names){
    stop("You must have a column containing the replicate measurement IDs, e.g., ReplicateID")
  }
  else if(!all(c("MP_A","MP_B") %in% data_names)){
    stop("You must have your data on LFDT format with columns MP_A and MP_B")
  }
  groups <- data_names[which(!data_names %in% c("MP_A","MP_B","ReplicateID") )]
  if(measure[1] == "standard deviation" | measure[1] == 1 | tolower(measure[1]) == "sd"){
    data <- data[, list(SD_A = sd(MP_A),
                        SD_B = sd(MP_B),
                        Concentration = (mean(MP_A) + mean(MP_B)) / 2), by = groups]
    data <- melt(data = data,
                 id.vars = c(groups,"Concentration"),
                 measure.vars = c("SD_A", "SD_B"),
                 variable.name = "Measure",
                 value.name = "Value",
                 variable.factor = FALSE)
    #### ok ####
    data <- data[, list(Measure = stri_split_fixed(str = Measure, pattern = "_", simplify = TRUE)[,1],
                        MS = stri_split_fixed(str = Measure, pattern = "_", simplify = TRUE)[,2],
                        Concentration = Concentration,
                        Value = Value), by = groups]
    data <- data[, list(Measure = Measure,
                        MS = ifelse(MS=="A",
                                    stri_split_fixed(Comparison, pattern = " - ", simplify = TRUE)[,1],
                                    stri_split_fixed(Comparison, pattern = " - ", simplify = TRUE)[,2]),
                        Concentration = Concentration,
                        Value = Value), by = groups]
    return(data)

  }
  else if(measure[1] == "coefficient of variation" | measure[1] == 2 | tolower(measure[1]) == "cv"){
    data <- data[, list(CV_A = sd(MP_A) / mean(MP_A),
                        CV_B = sd(MP_B) / mean(MP_B),
                        Concentration = (mean(MP_A) + mean(MP_B)) / 2), by = groups]
    data <- melt(data = data,
                 id.vars = c(groups,"Concentration"),
                 measure.vars = c("CV_A", "CV_B"),
                 variable.name = "Measure",
                 value.name = "Value",
                 variable.factor = FALSE)
    #### ok ####
    data <- data[, list(Measure = stri_split_fixed(str = Measure, pattern = "_", simplify = TRUE)[,1],
                        MS = stri_split_fixed(str = Measure, pattern = "_", simplify = TRUE)[,2],
                        Concentration = Concentration,
                        Value = Value), by = groups]
    data <- data[, list(Measure = Measure,
                        MS = ifelse(MS=="A",
                                    stri_split_fixed(Comparison, pattern = " - ", simplify = TRUE)[,1],
                                    stri_split_fixed(Comparison, pattern = " - ", simplify = TRUE)[,2]),
                        Concentration = Concentration,
                        Value = Value), by = groups]
    return(data)
  }
  else if(measure[1] == "lambda" | measure[1] == 3 | tolower(measure[1]) == "l"){
    data <- data[, list(lambda = var(MP_A) / var(MP_B),
                        Concentration = (mean(MP_A) + mean(MP_B)) / 2), by = groups]
    data <- melt(data = data,
                 id.vars = c(groups,"Concentration"),
                 measure.vars = c("lambda"),
                 variable.name = "Measure",
                 value.name = "Value",
                 variable.factor = FALSE)
    return(data)
  }
  else{
    message(paste0("measure = '",measure[1],"' is not a valid input. Please use method = 'standard deviation' or method = 'coefficient of variation'. SD is used as default..."))
    data <- data[, list(SD_A = sd(MP_A),
                        SD_B = sd(MP_B),
                        Concentration = (mean(MP_A) + mean(MP_B)) / 2), by = groups]
    data <- melt(data = data,
                 id.vars = c(groups,"Concentration"),
                 measure.vars = c("SD_A", "SD_B"),
                 variable.name = "Measure",
                 value.name = "Value",
                 variable.factor = FALSE)
    #### ok ####
    data <- data[, list(Measure = stri_split_fixed(str = Measure, pattern = "_", simplify = TRUE)[,1],
                        MS = stri_split_fixed(str = Measure, pattern = "_", simplify = TRUE)[,2],
                        Concentration = Concentration,
                        Value = Value), by = groups]
    data <- data[, list(Measure = Measure,
                        MS = ifelse(MS=="A",
                                    stri_split_fixed(Comparison, pattern = " - ", simplify = TRUE)[,1],
                                    stri_split_fixed(Comparison, pattern = " - ", simplify = TRUE)[,2]),
                        Concentration = Concentration,
                        Value = Value), by = groups]
    return(data)
  }
}

