#' Estimate prediction band using Deming regression formulated by J. Gillard and Fuller
#'
#' @param data A data frame or data table with format LFDT enclosed with all replicated measurements
#' @param level A numeric value that captures the overall confidence level of the estimated prediction band. This is assumed to already be Bonferroni-corrected. The recommended base level is 0.99
#' @param R An integer signifying the maximum number of replicated measurements performed on each evaluated material. This will shrink the width of the prediction interval and band because using the mean of replicates decreases uncertainty
#' @param Np An integer, which captures the number of pointwise prediction intervals making the prediction band across the concentration range. Ignored if evaluated_materials is !NULL
#' @param evaluated_materials A data frame or data table with format LFDT enclosed with all replicated measurements for evaluated materials such as EQAMs or CRMs. Should be NULL if the PB, that is, pointwise prediction intervals to be estimated
#'
#' @return A data table enclosed with information regarding the estimated prediction band across the concentration range or for the particular evaluated materials
#' @export
#'
#' @examples Deming_Gillard_Fuller(data = MS_wise(sampled_cs_measurements)[Comparison=="MP1 - MP2",], level = 0.95, R = 3, evaluated_materials = NULL)

Deming_Gillard_Fuller <- function(data, level = 0.99, R = 3, Np = 1e3, evaluated_materials = NULL){
  data <- as.data.table(data)
  R <- ceiling(R)
  Np <- ceiling(Np)
  precisions <- MS_precision(data)
  data <- mean_of_replicates(data)
  lambda <- precisions$lambda[1]
  x <- data$MP_B
  y <- data$MP_A
  n <- nrow(data)
  nx <- seq(from = min(x), to = max(x), length.out = Np)
  if(!is.null(evaluated_materials)){
    if(!any(class(evaluated_materials)[1] == c("data.table", "data.frame", "tibble", "list", "matrix", "array"))){
      message("For maintainer (evaluated_materials): evaluated_materials is not of correct class")
      stop("evaluated_materials is not data frame, data table or NULL. Please make sure that evaluated_materials's class is one of the data.frame, data.table or NULL. See ?Deming_Gillard_Fuller")
    }
    if(!any(class(evaluated_materials)[1] == c("list", "matrix", "array"))){
      warning(paste0(class(evaluated_materials)[1], " is the class of data. If not properly defined it may cause problems or unwanted results. I recommend you to use data.table, data.frame or tibble instead"))
      if(is.null(colnames(evaluated_materials))){
        message("For maintainer (evaluated_materials): column names of a matrix / array / list must not be NULL !")
        stop("Column names of a matrix / array / list must not be NULL! Make sure they are named by given standards of LFDT")
      }
    }
    evaluated_materials <- as.data.table(evaluated_materials)
    evaluated_materials <- mean_of_replicates(evaluated_materials)
    nx <- evaluated_materials$MP_B
  }
  sxx <- crossprod(x - mean(x))[,] / (n - 1)
  syy <- crossprod(y - mean(y))[,] / (n - 1)
  sxy <- crossprod(x - mean(x), y - mean(y))[,] / (n - 1)
  a <- syy - lambda * sxx
  b <- sqrt((syy - lambda * sxx)**2 + 4 * lambda * sxy**2)
  c <- 2 * sxy
  b1 <- (a + b) / c
  b0 <- mean(y) - b1 * mean(x)
  s.uu <- (syy + lambda * sxx - b)/(2*lambda)
  s.vv <- crossprod(y - mean(y) - b1 * (x - mean(x))) / (n - 2)
  var.b1 <-  (((b1**2) / (n * sxy**2))) * (sxx * syy - sxy**2)
  var.b0 <- s.vv / n + mean(x)**2 * var.b1
  cov.b0.b1 <- -mean(x) * var.b1
  v <- rbind(c(var.b0, cov.b0.b1), c(cov.b0.b1, var.b1))
  t <- qt(p = (1 - level) / 2, df = n - 2, lower.tail = FALSE)
  b <- c(b0, b1)
  X <- cbind(rep(1, length(nx)), nx)
  ny <- (X %*% b)[,]
  epsilon.model <- mean((lambda / (lambda + b1**2)) * x + (b1 / (lambda + b1**2)) * (y - b0))
  epsilon <- (lambda / (lambda + b1**2)) * nx + (b1 / (lambda + b1**2)) * (ny - b0)
  prediction.variance <- var.b1 * (epsilon - epsilon.model)**2 + (1 + 1 / n) * (s.uu / R) * (lambda + b1**2 + (var.b1 / (1 + 1 / n)))
  prediction.std.error <- sqrt(prediction.variance)
  return(data.table("X.GRID" = nx, "fit" = ny, "lwr" = ny - t * prediction.std.error, "upr" = ny + t * prediction.std.error))
}
