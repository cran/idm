\name{transform_z}
\alias{transform_z}
\title{
Computes indicator and standardized residual matrices
}
\description{
This function returns the indicator (or dummy) matrix and the standardized residual matrix of a given data matrix.
}
\usage{
transform_z(data,is.weight=TRUE,is.exact=TRUE,r=1,c=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Matrix or data frame}
  \item{is.weight}{Logical indicating whether computations are based on user-defined weights/masses}
  \item{is.exact}{Logical indicating whether the computations mimic the case when 'exact' or 'live' MCA is applied}
  \item{r}{Row masses; applicable only when is.weight = TRUE}
  \item{c}{Column masses; applicable only when is.weight = TRUE}
}
\value{
\item{dZ}{Indicator matrix}
\item{Q}{Number of variables}
\item{J}{Number of variable attributes}
\item{r}{Row masses}
\item{c}{Column masses}
\item{SZ}{Standardized residual matrix}
}
\examples{
data(women, package="idm")
outZ = transform_z(women[,c(1:7)],is.weight=FALSE)
}
