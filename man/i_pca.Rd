\name{i_pca}
\alias{i_pca}
\title{
Incremental Principal Component Analysis (PCA)
}
\description{
This function computes the Principal Component Analysis (PCA) solution on the covariance matrix using the incremental method of Hall, Marshall & Martin (2002).
}
\usage{
i_pca(data1, data2, nchunk = 2, disk = FALSE)
}
\arguments{
  \item{data1}{Matrix or data frame of starting data, or full data if data2 = NULL}
  \item{data2}{Matrix or data frame of incoming data; omitted when full data is given in data1}
  \item{nchunk}{Number of incoming data chunks (equal splits of 'data2') or a Vector with the row size of each incoming data chunk}
\item{disk}{Logical indicating whether then output is saved to hard disk}
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
\item{inertia_e}{Percentage of explained variance}
\item{sv}{Singular values}
\item{levelnames}{Column labels}
\item{rowctr}{Row contributions}
\item{colctr}{Column contributions}
\item{rowcor}{Row squared correlations}
\item{colcor}{Column squared correlations}
\item{nchunk}{Number of incoming data chunks}
\item{disk}{Logical indicating whether the output is saved to hard disk}
\item{allrowcoords}{A list containing the row scores on the principal components produced after each data chunk is analyzed; applicable only when disk = FALSE}
\item{allcolcoords}{A list containing the variable loadings on the principal components produced after each data chunk is analyzed; applicable only when disk = FALSE}
\item{allrowctr}{A list containing the row contributions after each data chunk is analyzed; applicable only when disk = FALSE}
\item{allrowcor}{A list containing the row squared correlations produced after each data chunk is analyzed; applicable only when disk = FALSE}
}
\references{
Hall, P., Marshall, D., & Martin, R. (2002). Adding and subtracting eigenspaces with eigenvalue decomposition and singular value decomposition. \emph{Image and Vision Computing}, \emph{20}(13), 1009-1016. \cr \cr
Iodice D' Enza, A., & Markos, A. (2015). Low-dimensional tracking of association structures in categorical data, \emph{Statistics and Computing}, \emph{25}(5), 1009-1022. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{update.i_pca}}, \code{\link{i_mca}}, \code{\link{update.i_mca}}, \code{\link{add_es}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#Iris Species 
data("iris", package = "datasets")
X = scale(iris[,-5])
res_iPCA = i_pca(data1 = X[1:50,], data2 = X[51:150,], nchunk = c(50,50))
#static plot
plot(res_iPCA, animation = FALSE) 
\donttest{
#\donttest is used here because the code calls the saveGIF function of the animation package 
#which requires ImageMagick or GraphicsMagick to be installed in your system 
#See help(im.convert) for details on the configuration of ImageMagick or GraphicsMagick.
#Creates animated GIF movies for objects and variables
plot(res_iPCA, animation = TRUE, frames = 10)
}

#Daily Closing Prices of Major European Stock Indices, 1991-1998 
data("EuStockMarkets", package = "datasets") 
res_iPCA = i_pca(data1 = EuStockMarkets[1:50,], data2 = EuStockMarkets[51:1860,], nchunk = 5) 
\donttest{
#\donttest is used here because the code calls the saveGIF function of the animation package 
#which requires ImageMagick or GraphicsMagick to be installed in your system 
#See help(im.convert) for details on the configuration of ImageMagick or GraphicsMagick.
#Creates animated GIF movies for objects and variables
plot(res_iPCA, animation = TRUE, frames = 10)
}
}