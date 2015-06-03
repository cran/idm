\name{update_pca}
\alias{update_pca}
\title{
Updates a Principal Component Analysis solution
}
\description{
This function updates the Principal Component Analysis (PCA) solution on the covariance matrix using the incremental method of Hall, Marshall & Martin (2002)
}
\usage{
update_pca(eg, incdata)
}

\arguments{
\item{incdata}{
Matrix of incoming data
}
\item{eg}{
A list describing the eigenspace of a data matrix, with components \cr
\code{vctCol} {Left eigenvectors} \cr
\code{vct} {Right eigenvectors} \cr
\code{n} {Number of cases} \cr
\code{val} {Eigenvalues} \cr
\code{orgn} {Data mean}
}
}

\value{
A list describing the eigenspace of a data matrix and a PCA solution, with components 
\item{vctCol}{Left eigenvectors} 
\item{vct}{Right eigenvectors} 
\item{n}{Number of cases} 
\item{val}{Eigenvalues}
\item{orgn}{Data mean}
\item{colpcoord}{Variable loadings}     
\item{rowpcoord}{Row scores on the principal components}     
\item{levelnames}{Column names}
\item{rowctr}{Row contributions}
\item{rowcor}{Row squared correlations}
\item{sv}{Eigenvalues}
\item{inertia_e}{Percentages of explained variance}
}
\references{
Hall, P., Marshall, D., & Martin, R. (2002). Adding and subtracting eigenspaces with eigenvalue decomposition and singular value decomposition. \emph{Image and Vision Computing}, 20(13), 1009-1016. \cr
Iodice D' Enza, A., & Markos, A. (2014). Low-dimensional tracking of association structures in categorical data. \emph{Statistics and Computing} (forthcoming).
}

\seealso{
\code{\link{update_mca}}, \code{\link{i_pca}}, \code{\link{i_mca}}, \code{\link{add_svd}}, \code{\link{add_eig}}
}

\examples{
data(iris, package="datasets")
dat = iris[,-5]
eg = do_eig(dat[1:50,-5])
aa = seq(50,150,50)
for (k in c(1:2)) {
  eg = update_pca(eg, dat[c((aa[k]+1):aa[k+1]),])
  }
plot(eg, animation = FALSE)  
}
