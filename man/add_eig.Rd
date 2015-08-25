\name{add_eig}
\alias{add_eig}
\title{
Adds two eigenspaces using incremental eigendecomposition
}
\description{
This function adds two eigenspaces using the incremental method of Hall, Marshall & Martin (2002). The results correspond to the batch SVD of the mean-centered and concatenated data.
}

\usage{
add_eig(m1, m2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{m1}{A list describing the eigenspace of a data matrix, with components \cr
\code{u} {Left eigenvectors} \cr
\code{v} {Right eigenvectors} \cr
\code{m} {Number of cases} \cr
\code{d} {Eigenvalues} \cr
\code{orgn} {Data mean} \cr
}
\item{m2}{A list describing the eigenspace of a data matrix, with components \cr
\code{u} {Left eigenvectors} \cr
\code{v} {Right eigenvectors} \cr
\code{m} {Number of cases} \cr
\code{d} {Eigenvalues} \cr
\code{orgn} {Data mean} \cr}
}

\value{
A list describing the eigenspace of a data matrix, with components 
\item{u}{Left eigenvectors} 
\item{v}{Right eigenvectors} 
\item{m}{Number of cases} 
\item{d}{Eigenvalues} 
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
X = iris[,-5]
#obtain two eigenspaces
m1 = do_eig(X[1:50,])
m2 = do_eig(X[c(51:150),])
#add the two eigenspaces keeping track of the data mean
m12 = add_eig(m1, m2)
#similar to the SVD of the mean-centered X
}

