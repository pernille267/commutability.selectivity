#' Estimate prediction band or prediction intervals using ordinary least squares regression
#'
#' @param data A data frame or data table with format LFDT enclosed with all replicated measurements
#' @param level A numeric value that captures the overall confidence level of the estimated prediction band. This is assumed to already be Bonferroni-corrected. The recommended base level is 0.99
#' @param R An integer signifying the maximum number of replicated measurements performed on each evaluated material. This will shrink the width of the prediction interval and band because using the mean of replicates decreases uncertainty
#' @param Np An integer, which captures the number of pointwise prediction intervals making the prediction band across the concentration range. Ignored if evaluated_materials is not \code{NULL}
#' @param evaluated_materials A data frame or data table with format LFDT enclosed with all replicated measurements for evaluated materials such as EQAMs or CRMs. Should be NULL if the PB, that is, pointwise prediction intervals to be estimated
#' @param optimize_for_measurement_errors Logical - Should we shift axis if the measurement error for MP_B is larger than MP_A? This will be optimal because OLS assumes that variability in MP_B is zero
#'
#' @return A data table enclosed with information regarding the estimated prediction band across the concentration range or for the particular evaluated materials
#' @export
#'
#' @examples Ordinary_Least_Squares(data = sdwinsp(R = 5)[,Comparison:="NO"], level = 0.90, R = 5, evaluated_materials = NULL)
Ordinary_Least_Squares <- function(data, level = 0.99, R = 3, Np = 1e3, evaluated_materials = NULL, optimize_for_measurement_errors = TRUE){
  data <- as.data.table(data)
  R <- ceiling(R)
  Np <- ceiling(Np)
  precisions <- MS_precision(data)
  data <- mean_of_replicates(data)
  lambda <- precisions$lambda[1]
  if(optimize_for_measurement_errors){
    if(lambda > 1){
      x <- data$MP_B
      y <- data$MP_A
      n <- nrow(data)
    }
    else{
      x <- data$MP_A
      y <- data$MP_B
      n <- nrow(data)
    }
  }
  else{
    x <- data$MP_B
    y <- data$MP_A
    n <- nrow(data)
  }
  nx <- seq(from = min(x), to = max(x), length.out = Np)
  if(!is.null(evaluated_materials)){
    evaluated_materials <- as.data.table(evaluated_materials)
    evaluated_materials <- mean_of_replicates(evaluated_materials)
    if(all(optimize_for_measurement_errors, lambda < 1)){
      nx <- evaluated_materials$MP_A
    }
    else{
      nx <- evaluated_materials$MP_B
    }

  }
  sxx <- crossprod(x - (sum(x)/length(x)))[,]/(n-1)
  syy <- crossprod(y - (sum(y)/length(y)))[,]/(n-1)
  sxy <- crossprod(x - (sum(x)/length(x)),y - (sum(y)/length(y)))[,]/(n-1)
  b1 <- sxy / sxx
  b0 <- sum(y) / length(y) - b1 * (sum(x) / length(x))
  s2 <- crossprod(y - (b0 + b1 * x))[,]/(n-2)
  t <- qt(p = (1-level)/2, df = n-2, lower.tail = FALSE)
  ny <- b0 + b1 * nx
  prediction_interval_width <- t * sqrt(s2 * (1 + 1/n + ((nx - mean(x))^2)/sxx) / R)
  return(data.table("X.GRID" = nx, "fit" = ny, "lwr" = ny - prediction_interval_width, "upr" = ny + prediction_interval_width, "reversed" = all(lambda < 1, optimize_for_measurement_errors)))
}
