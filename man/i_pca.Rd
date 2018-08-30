\name{i_pca}
\alias{i_pca}
\title{
Incremental Principal Component Analysis (PCA)
}
\description{
This function computes the Principal Component Analysis (PCA) solution on the covariance matrix using the incremental method of Hall, Marshall & Martin (2002).
}
\usage{
i_pca(data1, data2, current_rank, nchunk = 2, disk = FALSE)
}
\arguments{
  \item{data1}{Matrix or data frame of starting data, or full data if data2 = NULL}
  \item{data2}{Matrix or data frame of incoming data; omitted when full data is given in data1}
  \item{current_rank}{Rank of approximation or number of components to compute; if empty, the full rank is used}
  \item{nchunk}{Number of incoming data chunks (equal splits of 'data2', \code{default = 2}) or a Vector with the row size of each incoming data chunk}
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
\item{sv}{Singular values}
\item{inertia_e}{Percentage of explained variance}
\item{levelnames}{Attribute labels}
\item{rowctr}{Row contributions}
\item{colctr}{Column contributions}
\item{rowcor}{Row squared correlations}
\item{colcor}{Column squared correlations}
\item{nchunk}{A copy of \code{nchunk} in the return object}
\item{disk}{A copy of \code{disk} in the return object}
\item{allrowcoord}{A list containing the row scores on the principal components produced after each data chunk is analyzed; returned only when \code{disk = FALSE}}
\item{allcolcoord}{A list containing the variable loadings on the principal components produced after each data chunk is analyzed; returned only when \code{disk = FALSE}}
\item{allrowctr}{A list containing the row contributions after each data chunk is analyzed; returned only when \code{disk = FALSE}}
\item{allcolctr}{A list containing the column contributions after each data chunk is analyzed; returned only when \code{disk = FALSE}}
\item{allrowcor}{A list containing the row squared correlations produced after each data chunk is analyzed; returned only when \code{disk = FALSE}}
\item{allcolcor}{A list containing the column squared correlations produced after each data chunk is analyzed; returned only when \code{disk = FALSE}}
}
\references{
Hall, P., Marshall, D., & Martin, R. (2002). Adding and subtracting eigenspaces with eigenvalue decomposition and singular value decomposition. \emph{Image and Vision Computing}, \emph{20}(13), 1009-1016.\cr \cr
Iodice D' Enza, A., & Markos, A. (2015). Low-dimensional tracking of association structures in categorical data, \emph{Statistics and Computing}, \emph{25}(5), 1009--1022. \cr \cr
Iodice D'Enza, A., Markos, A., & Buttarazzi, D. (2018). The {idm} Package: Incremental Decomposition Methods in {R}. \emph{Journal of Statistical Software, Code Snippets}, \emph{86}(4), 1--24. DOI: 10.18637/jss.v086.c04.}

\seealso{
\code{\link{update.i_pca}}, \code{\link{i_mca}}, \code{\link{update.i_mca}}, \code{\link{add_es}}
}
\examples{
data("segmentationData", package = "caret")
#center and standardize variables, keep 58 continuous attributes
HCS = data.frame(scale(segmentationData[,-c(1:3)]))
#abbreviate variable names for plotting
names(HCS) = abbreviate(names(HCS), minlength = 5)
#split the data into starting data and incoming data
data1 = HCS[1:150, ]
data2 = HCS[151:2019, ]
#Incremental PCA on the HCS data set: the incoming data is 
#splitted into twenty chunks; the first 5 components/dimensions 
#are computed in each update
res_iPCA = i_pca(data1, data2, current_rank = 5, nchunk = 20)
#Static plots 
plot(res_iPCA, animation = FALSE)
\donttest{
#\donttest is used here because the code calls the saveLatex function of the animation package 
#which requires ImageMagick or GraphicsMagick and 
#Adobe Acrobat Reader to be installed in your system 
#See help(im.convert) for details on the configuration of ImageMagick or GraphicsMagick.
#Creates animated plot in PDF for objects and variables
plot(res_iPCA, animation = TRUE, frames = 10, movie_format = 'pdf')
}

#Daily Closing Prices of Major European Stock Indices, 1991-1998 
data("EuStockMarkets", package = "datasets") 
res_iPCA = i_pca(data1 = EuStockMarkets[1:50,], data2 = EuStockMarkets[51:1860,], nchunk = 5) 
\donttest{
#\donttest is used here because the code calls the saveLatex function of the animation package 
#which requires ImageMagick or GraphicsMagick and 
#Adobe Acrobat Reader to be installed in your system 
#See help(im.convert) for details on the configuration of ImageMagick or GraphicsMagick.
#Creates animated plot in PDF movies for objects and variables
plot(res_iPCA, animation = TRUE, frames = 10, movie_format = 'pdf')
}
}