#' Studentize data using the standard technique
#'
#' @param data A numeric vector that we want to studentize
#'
#' @return A numeric vector that has studentized values
#' @export
#'
#' @examples studentize(rnorm(100,1,1))

studentize <- function(data){
  return((data - mean(data, na.rm = TRUE))/sd(data, na.rm = TRUE))
}
