% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MS_wise_2.R
\name{MS_wise_2}
\alias{MS_wise_2}
\title{Converting standard wide format MS comparison data (WFDT) to standard long format (LFDT)}
\usage{
MS_wise_2(data, fix_names = TRUE)
}
\arguments{
\item{data}{A data frame or data table with m + n columns where m (usually two) of them are ID columns and n of them is numeric columns consisting of measurement results for each measurement system. The ID columns may be named approximately e.g., SampleID or ReplicateID, but names deviating much from these will not be recongnized}

\item{fix_names}{Should ID column names with inappropriate names be attempted fixed. Default is \code{TRUE}, if we set \code{FALSE} and the ID columns not named correctly, this function will produce an error}
}
\value{
A data table of format LFDT, with m + 3 columns where the latter three columns are MS comparison ID column (new), Sample ID column and Replicate ID column. The m remaining columns are the resulting numeric columns representing MS measurements
}
\description{
A backup function or convenience function where we do not need to specify the ID columns explicitly. The given ID column names must however be close to corresponding standard names
}
\examples{
MS_wise_2(sampled_cs_measurements, fix_names = FALSE)
}
