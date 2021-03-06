% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/present_MS_precision_estimates.R
\name{present_MS_precision_estimates}
\alias{present_MS_precision_estimates}
\title{Present precision estimates for unique measurement systems}
\usage{
present_MS_precision_estimates(
  data,
  groups = "Comparison",
  LFDT = FALSE,
  decimals = 6
)
}
\arguments{
\item{data}{A data frame or data table with format LFDT or WFDT enclosed with all replicated measurements}

\item{groups}{Character vector - The names of the grouping columns of the data. Default is 'Comparison'}

\item{LFDT}{Logical - Is data of format LDFT or WFDT (default)?}

\item{decimals}{Integer - how many decimals should be used in the output data table. Default is 6}
}
\value{
A data table with two rows and 1 + n columns, where n is the number of unique measurement systems in data
}
\description{
Present precision estimates for unique measurement systems
}
\examples{
present_MS_precision_estimates(sampled_cs_measurements)
}
