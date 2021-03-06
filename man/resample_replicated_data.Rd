% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resample_replicated_data.R
\name{resample_replicated_data}
\alias{resample_replicated_data}
\title{Resample data with replicated measurements}
\usage{
resample_replicated_data(data)
}
\arguments{
\item{data}{A data frame or data table with format LFDT where all replicate measurements are enclosed. Two of the columns must be the standard ID columns SampleID and ReplicateID}
}
\value{
A data table that is the resulting resampled data with CSs as clusters
}
\description{
Resample data with replicated measurements
}
\examples{
resample_replicated_data(data = MS_wise(sampled_cs_measurements)[Comparison=="MP1 - MP3",])
}
