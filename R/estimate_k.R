#' Estimate a k value based on one MS comparison
#'
#' @param data A data frame or data table with format LFDT where all replicate measurements are enclosed
#'
#' @return A data table containing a plug-in point estimate of K, namely, k
#' @export
#'
#' @examples estimate_k(data = MS_wise(sampled_cs_measurements)[Comparison=="MP1 - MP2"])

estimate_k <- function(data){
  names <- names(data)[!names(data)%in%c("ReplicateID","MP_A","MP_B")]
  N <- length(data$MP_A)
  mpprec <- MS_precision(data = data)
  if(mpprec$lambda<1){
    mod <- lm.fit(x=cbind(rep(1,N),data$MP_A),y=data$MP_B)
    var.res <- crossprod(x=mod$residuals,y=mod$residuals)[,]/mod$df.residual
    b1 <- mod$coefficients[2]
    bias <- (N+2)/N
    k <- (bias*var.res)/(mpprec$Var_B + mpprec$Var_A * (b1^2))
    return(data.table(k=k))
  }
  mod <- lm.fit(x=cbind(rep(1,N),data$MP_A),y=data$MP_B)
  var.res <- crossprod(x=mod$residuals,y=mod$residuals)[,]/mod$df.residual
  b1 <- mod$coefficients[2]
  bias <- (N+2)/N
  k <- (bias*var.res)/(mpprec$Var_B * (b1^2) + mpprec$Var_A)
  return(data.table(k=k))
}



