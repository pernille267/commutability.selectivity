#' Based on lambda, which set of measurement system measurements should be used as x and y
#'
#' @param data A data table found in \code{Ordinary_Least_Squares_over_groups()}
#' @param groups A vector of grouping columns
#'
#' @return A fixed data table
#' @export
#'
#' @examples fix_reverse_issues(data.table(Comparison = "f - g", SampleID = "ID 1", MP_B = 5.10, MP_A = 5.20, fit = 5.08, lwr = 5.00, upr = 5.16, reversed = TRUE))
fix_reverse_issues <- function(data, groups = "Comparison"){
  data[,.(`Comparison order` = ifelse(reversed, paste0(rev(unlist(stri_split(str = Comparison, fixed = " - "))),collapse=" - "), Comparison),
          SampleID = SampleID,
          MP_B = MP_B,
          MP_A = MP_A,
          x = ifelse(reversed,MP_A,MP_B),
          y = ifelse(reversed,MP_B,MP_A),
          fit = fit,
          lwr = lwr,
          upr = upr), by = groups]
}

