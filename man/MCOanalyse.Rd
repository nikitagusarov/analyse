% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mco_selection.r
\name{MCOanalyse}
\alias{MCOanalyse}
\title{MCO automation function : Parallelized approach}
\usage{
MCOanalyse(datax, dependent, laglim = 2, plim = 0.1, limit = 1,
  study.t0 = FALSE, varlim = 10, valid = TRUE, part = 0.15)
}
\arguments{
\item{datax}{A data.frame of explicative variables to be studied of length equal to explicative variable vector}

\item{dependent}{A dependent variable vector}

\item{laglim}{The model allows to automatically convert the dataset in order to study lagged dependecies}

\item{plim}{The critical value for tests results (default value 0.1)}

\item{limit}{A maximum number of models to be selected in each category}

\item{study.t0}{A logical value indicating whether the t0 period should be taken into account}

\item{varlim}{A maximal number of variables to make combinations of (large number significantly increases computation time, each supplementary variable increases the tyme by 2)}

\item{valid}{LOGICAL, indicates whether the crossvalidation procedure should be applied during the model selection}

\item{part}{The part of the sample to be used for crossvalidation (by default 15%)}
}
\value{
A list of selected models estimated with OLS with additional performance information
\describe{
 \item{choices}{A data.frame containing information for selected models testing results}
 \item{models}{A list of models in order as they figure in choices dataframe}
 \item{univariate}{Results of first stage univariate model testing for selected variables}
}
}
\description{
This function allows to select by informational criteria and various tests results (Durbin-Watson, Shapiro-Wilks, Breush-Pagan, F-test) one or a serie of different models allowing for best performance in-sample. The function uses parallelised approach with doParallel package.
}
\examples{
MCOanalyse(data, dependent, 
     laglim = 5, plim = 0.1, limit = 2, 
     study.t0 = TRUE, folds = 10, varlim = 10)
}
\keyword{MCOanalyse}
