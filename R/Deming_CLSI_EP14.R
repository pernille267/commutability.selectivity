#' Estimate prediction band using Deming regression formulated as in EP14
#'
#' @param data A data frame or data table with format LFDT enclosed with all replicated measurements
#' @param level A numeric value that captures the overall confidence level of the estimated prediction band. This is assumed to be NOT Bonferroni corrected. Recommended level is 0.99
#' @param R An integer signifying the maximum number of replicated measurements performed on each evaluated material. This will shrink the width of the prediction interval and band because using the mean of replicates decreases uncertainty
#' @param Np An integer, which captures the number of pointwise prediction intervals making the prediction band across the concentration range. Ignored if evaluated_materials is !NULL
#' @param evaluated_materials A data frame or data table with format LFDT enclosed with all replicated measurements for evaluated materials such as EQAMs or CRMs. Should be NULL if the PB, that is, pointwise prediction intervals to be estimated
#'
#' @return A data table enclosed with information regarding the estimated prediction band across the concentration range
#' @export
#'
#' @examples Deming_CLSI_EP14(MS_wise(sampled_cs_measurements)[Comparison=="MP1 - MP2",], level = 0.95, R = 3, Np = 3e3, evaluated_materials = NULL)

Deming_CLSI_EP14 <- function(data, level = 0.99, R = 3, Np = 1e3, evaluated_materials = NULL){
  R <- ceiling(R)
  Np <- ceiling(Np)
  precision <- MS_precision(data)
  data <- mean_of_replicates(data = data)
  lambda <- precision$lambda
  var.mp1 <- precision$Var_B
  x <- data$MP_B; y <- data$MP_A; n <- nrow(data)
  nx <- seq(from = min(x), to = max(x), length.out = Np)
  if(!is.null(evaluated_materials)){
    evaluated_materials <- as.data.table(evaluated_materials)
    evaluated_materials <- mean_of_replicates(evaluated_materials)
    nx <- evaluated_materials$MP_B
  }
  sxx <- (1/(n-1))*crossprod(x - mean(x))[,]
  syy <- (1/(n-1))*crossprod(y - mean(y))[,]
  sxy <- (1/(n-1))*crossprod(x - mean(x), y - mean(y))[,]
  b1 <- ((syy - lambda * sxx) + sqrt((syy - lambda * sxx)^2 + 4*lambda*sxy^2))/(2*sxy)
  b0 <- mean(y) - b1 * mean(x)
  s.uu <- var.mp1
  s.ee <- lambda * s.uu
  var.b1 <- ((b1 ^ 2 / (n * sxy^2))) * (sxx*syy - sxy^2)
  s.pred <- sqrt(((nx - mean(x))^2) * var.b1 + (1+1/n) * ((b1^2) * s.uu + s.ee)/R)
  t <- qt(p = (1-level)/2, df = n*(R-1), lower.tail = FALSE)
  b.vector <- c(b0,b1)
  x.new <- cbind(rep(1,times=length(nx)),nx)
  ny <- (x.new %*% b.vector)[,]
  predictions.clsi <- data.table("X.GRID" = nx, "fit" = ny, "lwr" = ny - t*s.pred, "upr" = ny + t*s.pred)
  return(predictions.clsi)
}
