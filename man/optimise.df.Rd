% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimise_df.r
\name{optimise.df}
\alias{optimise.df}
\title{Seasonal data.frames optimisation function}
\usage{
optimise.df(datax, season)
}
\arguments{
\item{datax}{A data.frame of seasonal data to be studied of length equal to explicative variable vector. The first column is advised to be year and the second to be "M" for monthly data, or "T" for trimestrial data.}

\item{season}{An indicator of seasonality type to be used as present in the data.frame.}
}
\value{
A data.frame with test results by variable.
}
\description{
This function allows to optimise the time series seasonal data, based on the stationarity tests' results.
It is an extention of test.df() function.
}
\examples{
optimise.df = function(datax, season)

}
\keyword{optimise.df}
\keyword{test.df}
