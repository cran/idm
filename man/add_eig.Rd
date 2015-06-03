\name{add_eig}
\alias{add_eig}
\title{
Adds two eigenspaces using incremental eigendecomposition
}
\description{
This function adds two eigenspaces using the incremental method of Hall, Marshall & Martin (2002)
}

\usage{
add_eig(m1, m2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{m1}{A list describing the eigenspace of a data matrix, with components \cr
\code{vctCol} {Left eigenvectors} \cr
\code{vct} {Right eigenvectors} \cr
\code{n} {Number of cases} \cr
\code{val} {Eigenvalues} \cr
\code{orgn} {Data mean} \cr
}
\item{m2}{A list describing the eigenspace of a data matrix, with components \cr
\code{vctCol} {Left eigenvectors} \cr
\code{vct} {Right eigenvectors} \cr
\code{n} {Number of cases} \cr
\code{val} {Eigenvalues} \cr
\code{orgn} {Data mean} \cr}
}

\value{
A list describing the eigenspace of a data matrix, with components 
\item{vctCol}{Left eigenvectors} 
\item{vct}{Right eigenvectors} 
\item{n}{Number of cases} 
\item{val}{Eigenvalues} 
\item{orgn}{Data mean} 
}


\references{
Hall, P., Marshall, D., & Martin, R. (2002). Adding and subtracting eigenspaces with eigenvalue decomposition and singular value decomposition. \emph{Image and Vision Computing}, 20(13), 1009-1016.
}

\seealso{
\code{\link{add_svd}}, \code{\link{update_pca}}, \code{\link{i_pca}}, \code{\link{i_mca}}, \code{\link{update_mca}}, \code{\link{do_eig}}
}
\examples{
#Iris species
data("iris", package = "datasets")
#obtain two eigenspaces
m1 = do_eig(iris[c(1:50),-5])
m2 = do_eig(iris[c(51:150),-5])
#add eigenspaces
m12 = add_eig(m1, m2)
}

