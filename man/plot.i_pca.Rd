\name{plot.i_pca}
\alias{plot.i_pca}
\title{
Plotting 2D maps in Principal Component Analysis
}
\description{
Graphical display of Principal Component Analysis results in two dimensions
}
\usage{\method{plot}{i_pca}(x, dims = c(1,2), what = c(TRUE,TRUE), 
dataname = NULL, labels = NULL, animation = TRUE, frames = 10, zoom = TRUE, ...)}

\arguments{
  \item{x}{Principal component analysis object returned by \code{\link{i_pca}}}
    \item{dims}{Numerical vector of length 2 indicating the dimensions to plot on horizontal and vertical axes respectively; default is first dimension horizontal and second dimension vertical}
\item{what}{Vector of two logicals specifying the contents of the plot(s). First entry indicates if the rows (objects) are displayed and the second entry if the columns (attributes) are displayed}
\item{dataname}{String prefix used for custom naming of output files; default is the name of the output object}
  \item{labels}{String vector of variable labels}
  \item{animation}{Logical indicating whether animated GIF files are created and saved to the hard drive or a static plot is obtained}
\item{frames}{Number of animation frames shown per iteration; applicable only when animation = TRUE}
\item{zoom}{Logical indicating whether axis limits change during the animation creating a zooming effect; applicable only when animation = TRUE}
  \item{...}{Further arguments passed to \code{\link{plot}} and \code{\link{points}}}
}

\details{
The function \code{plot.i_pca} makes a two-dimensional map of the object created by \code{i_pca} with respect to two selected dimensions.  
}
\references{ImageMagick: \url{http://www.imagemagick.org}; GraphicsMagick: \url{http:
//www.graphicsmagick.org}
}

\seealso{
\code{\link{plot.i_mca}}
}
\examples{
data("iris", package = "datasets")
#Iris species
X = scale(iris[,-5])
res_iPCA = i_pca(data1 = X[1:50,-5], data2 = X[51:150,-5], nchunk = c(50,50))
#static plot, final solution
plot(res_iPCA, animation = FALSE)
\donttest{
##\donttest is used here because the code calls the saveGIF function of the animation package 
#which requires ImageMagick or GraphicsMagick to be installed in your system 
#See help(im.convert) for details on the configuration of ImageMagick or GraphicsMagick.
#Creates animated GIF movies for objects and variables
plot(res_iPCA, animation = TRUE, frames = 10)
}
}

