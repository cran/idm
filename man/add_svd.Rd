\name{add_svd}
\alias{add_svd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adds two eigenspaces using incremental SVD
}
\description{
This function adds two eigenspaces using the incremental method of Ross et al. (2008)
}
\usage{
add_svd(U,S,V,B,m,orgn,f,current_rank)
}

\arguments{
\item{U}{Left singular vectors}
\item{S}{Singular values}
\item{V}{Right singular vectors}
\item{B}{Data matrix}
\item{m}{Number of cases processed up to this point}
\item{orgn}{Data mean}
\item{f}{Number between 0 and 1 indicating the forgetting factor}
\item{current_rank}{Parameter used to maintain the rank of the solution}
}

\value{
A list describing the eigenspace of a data matrix, with components \cr
\item{u}{Left singular vectors}
\item{d}{Singular values}
\item{v}{Right singular vectors}
\item{m}{Number of cases}
\item{orgn}{Data mean}
}
\references{
Ross, D. A., Lim, J., Lin, R. S., & Yang, M. H. (2008). Incremental learning for robust visual tracking. \emph{International Journal of Computer Vision}, 77(1-3), 125-141.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{add_eig}}, \code{\link{update_mca}}, \code{\link{i_pca}}, \code{\link{i_mca}}
}
\examples{
#Iris species
data("iris", package = "datasets")
#obtain eigenspace of the first 50 Iris species
eg = do_eig(data.matrix(iris[c(1:50),-5]))
#add this eigenspace to the one of the remaining cases
eg = add_svd(eg$vctCol,eg$val,eg$vct,data.matrix(iris[c(51:150),-5]),50,eg$orgn,f=1,current_rank=4)
}

