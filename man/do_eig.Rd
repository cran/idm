\name{do_eig}
\alias{do_eig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes the eigenspace of a data matrix
}
\description{
This function computes the eigenspace of a data matrix
}
\usage{
do_eig(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a matrix or data frame}
}

\value{
A list describing the eigenspace of a data matrix, with components 
\item{vctCol}{Left eigenvectors}
\item{vct}{Right eigenvectors}
\item{n}{Number of cases}
\item{val}{Eigenvalues}
\item{orgn}{Data mean}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{add_svd}}, \code{\link{update_pca}}, \code{\link{i_pca}}, \code{\link{i_mca}}, \code{\link{do_eig}}
}
\examples{
#Iris species
data("iris", package = "datasets")
m = do_eig(iris[,-5])
}

