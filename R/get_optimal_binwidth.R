#' Calculate the optimal binwidths for histograms based on several methods
#'
#' @param data A numeric vector that is used to construct histograms
#' @param method Which method should be used to get the optimal binwidths. Valid inputs are \code{FDR}, \code{SNRR}, \code{SF}, \code{RR} and \code{DF}. See details for more information regarding the different approaches
#'
#' @return A value that represents the optimal binwidth based on data and chosen method
#' @export
#'
#' @examples get_optimal_binwidth(rnorm(100), method = "SF")

get_optimal_binwidth <- function(data, method = c("FDR", "SNRR", "SF", "RR", "DF")){
  n <- length(data)
  propval <- n**(1/3)
  if(method[1] == "FDR"){
    h <- 2 * IQR(x = data, na.rm = TRUE) / propval
    return(h)
  }
  else if(method[1] == "SNRR"){
    h <- 3.49 * sd(data) / propval
    return(h)
  }
  else if(method[1] == "SF"){
    k <- log2(x = n)
    h <- (max(data) - min(data)) / k
    return(h)
  }
  else if(method[1] == "RR"){
    k <- 2 * propval
    h <- (max(data) - min(data)) / k
    return(h)
  }
  else if(method[1] == "DF"){
    numerator <- sum((na.omit(data) - mean(data, na.rm = TRUE))**3) / n
    denominator <- sd(data, na.rm = TRUE) ** 3
    g1 <- abs(numerator / denominator)
    sigma_g1 <- (6 * (n - 2) / (n + 1) / (n + 3)) ** (1 / 2)
    k <- 1 + log2(n) + log2(1 + g1 / sigma_g1)
    h <- (max(data) - min(data)) / k
    return(h)
  }
  else{
    message("method is not recongnized. Make sure you have spelled it correctly. FDR used by default")
    h <- 2 * IQR(x = data, na.rm = TRUE) / propval
    return(h)
  }
}
