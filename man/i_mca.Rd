\name{i_mca}
\alias{i_mca}
\title{
Incremental Multiple Correspondence Analysis (MCA)
}
\description{
This function computes the Multiple Correspondence Analysis (MCA) solution on the indicator matrix using two incremental methods described in Iodice D'Enza & Markos (2015)
}
\usage{
i_mca(data1, data2, method=c("exact","live"), current_rank, nchunk = 2, 
 ff = 0, disk = FALSE)
}
\arguments{
  \item{data1}{Matrix or data frame of starting data or full data if data2 = NULL}
  \item{data2}{Matrix or data frame of incoming data}
  \item{method}{String specifying the type of implementation: \code{"exact"} or \code{"live"}. \code{"exact"} refers to the case when all the data is available from the start and dimension reduction is based on the method of Hall et al. (2002). \code{"live"} refers to the case when new data comes in as data flows and dimension reduction is based on the method of Ross et al. (2008). The main difference between the two approaches lies in the calculation of the column margins of the input matrix. For the \code{"exact"} approach, the analysis is based on the "global" margins, that is, the margins of the whole indicator matrix, which is available in advance. For the \code{"live"} approach, the whole matrix is unknown and the global margins are approximated by  the "local" margins, that is, the average margins of the data analysed insofar. A detailed description of the two implementations is provided in Iodice D' Enza & Markos (2015).}
  \item{current_rank}{Rank of approximation or number of components to compute; if empty, the full rank is used}
  \item{nchunk}{Number of incoming data chunks (equal splits of 'data2', \code{default = 2}) or a Vector with the row size of each incoming data chunk}
  \item{ff}{Number between 0 and 1 indicating the "forgetting factor" used to down-weight the contribution of earlier data blocks to the current solution. When \code{ff = 0} (default) no forgetting occurs; applicable only when \code{method ="live"}}
  \item{disk}{Logical indicating whether then output is saved to hard disk}
}

\value{
\item{rowpcoord}{Row principal coordinates}     
\item{colpcoord}{Column principal coordinates}
\item{rowcoord}{Row standard coordinates}     
\item{colcoord}{Column standard coordinates}
\item{sv}{Singular values}
\item{inertia.e}{Percentages of explained inertia}
\item{levelnames}{Column labels}
\item{rowctr}{Row contributions}
\item{colctr}{Column contributions} 
\item{rowcor}{Row squared correlations}
\item{colcor}{Column squared correlations}
\item{rowmass}{Row masses}
\item{colmass}{Column masses}  
\item{nchunk}{A copy of \code{nchunk} in the return object}
\item{disk}{A copy of \code{disk} in the return object}
\item{ff}{A copy of \code{ff} in the return object}
\item{allrowcoord}{A list containing the row principal coordinates  produced after each data chunk is analyzed; returned only when \code{disk = FALSE}}
\item{allcolcoord}{A list containing the column principal coordinates  on the principal components produced after each data chunk is analyzed; returned only when \code{disk = FALSE}}
\item{allrowctr}{A list containing the row contributions after each data chunk is analyzed; returned only when \code{disk = FALSE}}
\item{allcolctr}{A list containing the column contributions after each data chunk is analyzed; returned only when \code{disk = FALSE}}
\item{allrowcor}{A list containing the row squared correlations produced after each data chunk is analyzed; returned only when \code{disk = FALSE}}
\item{allcolcor}{A list containing the column squared correlations produced after each data chunk is analyzed; returned only when \code{disk = FALSE}}

}
\references{
Hall, P., Marshall, D., & Martin, R. (2002). Adding and subtracting eigenspaces with eigenvalue decomposition and singular value decomposition. \emph{Image and Vision Computing}, \emph{20}(13), 1009-1016.\cr \cr
Ross, D. A., Lim, J., Lin, R. S., & Yang, M. H. (2008). Incremental learning for robust visual tracking. \emph{International Journal of Computer Vision}, \emph{77}(1-3), 125-141.\cr \cr
Iodice D' Enza, A., & Markos, A. (2015). Low-dimensional tracking of association structures in categorical data, \emph{Statistics and Computing}, \emph{25}(5), 1009-1022. 
}

\seealso{
\code{\link{update.i_mca}}, \code{\link{i_pca}}, \code{\link{update.i_pca}}, \code{\link{add_es}}}

\examples{
##Example 1 - Exact case
data("women", package = "idm")
nc = 5 # number of chunks
res_iMCAh = i_mca(data1 = women[1:300,1:7], data2 = women[301:2107,1:7]
,method = "exact", nchunk = nc)
#static MCA plot of attributes on axes 2 and 3
plot(x = res_iMCAh, dim = c(2,3), what = c(FALSE,TRUE), animation = FALSE)
\donttest{
#\donttest is used here because the code calls the saveLatex function of the animation package 
#which requires ImageMagick or GraphicsMagick and 
#Adobe Acrobat Reader to be installed in your system 
#Creates animated plot in PDF for objects and variables
plot(res_iMCAh, animation = TRUE, frames = 10, movie_format = 'pdf')
}

##Example 2 - Live case
data("tweet", package = "idm")
nc = 5
#provide attributes with custom labels
labels = c("HLTN", "ICN", "MRT","BWN","SWD","HYT","CH", "-", "-/+", "+", "++", "Low", "Med","High")
#mimics the 'live' MCA implementation 
res_iMCAl = i_mca(data1 = tweet[1:100,], data2 = tweet[101:1000,],
method="live", nchunk = nc, current_rank = 2)

\donttest{
#\donttest is used here because the code calls the saveLatex function of the animation package 
#which requires ImageMagick or GraphicsMagick and 
#Adobe Acrobat Reader to be installed in your system 
#See help(im.convert) for details on the configuration of ImageMagick or GraphicsMagick.
#Creates animated plot in PDF for observations and variables
plot(res_iMCAl, labels = labels, animation = TRUE, frames = 10, movie_format = 'pdf')
}
}



