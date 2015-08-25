\name{plot.i_mca}
\alias{plot.i_mca}
\title{
Plotting 2D maps in Multiple Correspondence Analysis
}
\description{
Graphical display of Multiple Correspondence Analysis results in two dimensions
}
\usage{\method{plot}{i_mca}(x, dims=c(1,2), what=c(TRUE,TRUE), 
contrib="none", dataname=NULL, labels=NULL, animation=TRUE, frames=10, 
zoom=TRUE,...)
}

\arguments{
\item{x}{Multiple correspondence analysis object returned by \code{\link{i_mca}}}
\item{dims}{Numerical vector of length 2 indicating the dimensions to plot on horizontal and vertical axes respectively; default is first dimension horizontal and second dimension vertical}
\item{what}{Vector of two logicals specifying the contents of the plot(s). First entry indicates if the rows (objects) are displayed and the second entry if the columns (attributes) are displayed)}
\item{contrib}{Vector of two character strings specifying if attribute contributions should be represented by different label size. Available options are\cr
                 \kbd{"none"} (contributions are not indicated in the plot)\cr
                 \kbd{"cor"} (relative contributions are indicated by label size)\cr
                 \kbd{"ctr"} (absolute contributions are indicated by label size)\cr
The higher the contribution of a point, the larger its label size. Default is \kbd{"none"}}
\item{dataname}{String prefix used for custom naming of output files; default is the name of the output object}
\item{labels}{String vector of variable labels}
\item{animation}{Logical indicating whether animated GIF files are created and saved to the hard drive or a static plot is obtained}
 \item{frames}{Number of animation frames shown per iteration; applicable only when animation = TRUE}
 \item{zoom}{Logical indicating whether axis limits change during the animation creating a zooming effect; applicable only when animation = TRUE}
\item{...}{Further arguments passed to \code{\link{plot}} and \code{\link{points}}}
}
\details{
The function \code{plot.i_mca} makes a two-dimensional map of the object created by \code{i_mca} with respect to two selected dimensions. In this map both the row and column points are scaled to have inertias (weighted variances) equal to the principal inertia (eigenvalue or squared singular value) along the principal axes, that is both rows and columns are in pricipal coordinates. 
}

\references{
Greenacre, M.J. (1993) \emph{Correspondence Analysis in Practice}. London: Academic Press. \cr
Greenacre, M.J. (1993) Biplots in Correspondence Analysis, \emph{Journal of Applied Statistics}, 20, 251 - 269.\cr
ImageMagick: \url{http://www.imagemagick.org}; GraphicsMagick: \url{http:
//www.graphicsmagick.org}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot.i_pca}}
}
\examples{
data("women", package="idm")
nc = 4 # number of chunks
### Different approaches to incremental MCA:
# Live case
res_iMCAl = i_mca(data1 = women[1:50,1:4],data2 = women[51:300,1:4],method = "live",nchunk = nc)
#static plot, final solution
plot(res_iMCAl, contrib="ctr", animation = FALSE)
\donttest{
#\donttest is used here because the code calls the saveGIF function of the animation package 
#which requires ImageMagick or GraphicsMagick to be installed in your system 
#See help(im.convert) for details on the configuration of ImageMagick or GraphicsMagick.
#Creates animated GIF movies for objects and variables
plot(res_iMCAl, contrib = "ctr", animation = TRUE, frames = 10)
}
}
