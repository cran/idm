\name{do_es}
\alias{do_es}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes the eigenspace of a data matrix
}
\description{
This function computes the eigenspace of a mean-centered data matrix
}
\usage{
do_es(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a matrix or data frame}
}

\value{
A list describing the eigenspace of a data matrix, with components 
\item{u}{Left eigenvectors}
\item{v}{Right eigenvectors}
\item{m}{Number of cases}
\item{d}{Eigenvalues}
\item{orgn}{Data mean}
\item{smfq}{...}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{add_es}}, \code{\link{update.i_pca}}, \code{\link{i_pca}}
}

\examples{
#Iris species
data("iris", package = "datasets")
eg = do_es(iris[,-5]) 
#corresponds to the SVD of the centered data matrix
}

