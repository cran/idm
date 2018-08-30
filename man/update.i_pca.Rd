\name{update.i_pca}
\alias{update.i_pca}
\title{
Updates a Principal Component Analysis solution
}
\description{
This function updates the Principal Component Analysis (PCA) solution on the covariance matrix using the incremental method of Hall, Marshall & Martin (2002)
}
\usage{
\method{update}{i_pca}(object, incdata, current_rank, ...)
}

\arguments{
\item{object}{object of class 'i_pca'}
\item{incdata}{matrix of incoming data}
\item{current_rank}{Rank of approximation or number of components to compute; if empty, the full rank is used}
\item{...}{Further arguments passed to \code{\link{update}}}
}

\value{
\item{rowpcoord}{Row scores on the principal components}     
\item{colpcoord}{Variable loadings}     
\item{eg}{A list describing the eigenspace of a data matrix, with components \cr
\code{u} {Left eigenvectors} \cr
\code{v} {Right eigenvectors} \cr
\code{m} {Number of cases} \cr
\code{d} {Eigenvalues} \cr
\code{orgn} {Data mean} \cr}
\item{inertia.e}{Percentages of explained variance}
\item{sv}{Singular values}
\item{levelnames}{Variable names}
\item{rowcor}{Row squared correlations}
\item{rowctr}{Row contributions}
\item{colcor}{Column squared correlations}
\item{colctr}{Column contributions}
}
\references{
Hall, P., Marshall, D., & Martin, R. (2002). Adding and subtracting eigenspaces with eigenvalue decomposition and singular value decomposition. \emph{Image and Vision Computing}, \emph{20}(13), 1009-1016.\cr \cr
Iodice D' Enza, A., & Markos, A. (2015). Low-dimensional tracking of association structures in categorical data, \emph{Statistics and Computing}, \emph{25}(5), 1009--1022. \cr \cr
Iodice D'Enza, A., Markos, A., & Buttarazzi, D. (2018). The {idm} Package: Incremental Decomposition Methods in {R}. \emph{Journal of Statistical Software, Code Snippets}, \emph{86}(4), 1--24. DOI: 10.18637/jss.v086.c04.}

\seealso{\code{\link{update.i_mca}}, \code{\link{i_pca}}, \code{\link{i_mca}}, \code{\link{add_es}}}

\examples{
data(segmentationData, package = "caret")
HCS = data.frame(scale(segmentationData[,-c(1:3)]))
names(HCS) = abbreviate(names(HCS), minlength = 5)
res_PCA = i_pca(HCS[1:200, ])
aa = seq(from = 201, to = nrow(HCS), by = 200)
aa[length(aa)] = nrow(HCS)+1
for (k in c(1:(length(aa)-1))){
     res_PCA = update(res_PCA, HCS[c((aa[k]):(aa[k+1]-1)),])
    }
#Static plot
plot(res_PCA, animation = FALSE)
}
