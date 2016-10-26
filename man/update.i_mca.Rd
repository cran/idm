\name{update.i_mca}
\alias{update.i_mca}
\title{
Updates a Multiple Correspondence Analysis solution
}
\description{
This function updates the Multiple Correspondence Analysis (MCA) solution on the indicator matrix using the incremental method of Ross, Lim, Lin, & Yang (2008)
}
\usage{
\method{update}{i_mca}(object, incdata, current_rank, ff = 0, nchunk = 1, ...)
}

\arguments{
\item{object}{object of class 'i_mca'}
\item{incdata}{Matrix of incoming data}
\item{current_rank}{Rank of approximation; if empty then full rank is used}
\item{ff}{Number between 0 and 1 indicating the "forgetting factor" used to down-weight the contribution of earlier data blocks to the current solution. When \code{ff = 0} (default) no forgetting occurs}
\item{nchunk}{Number of data blocks processed up to that point; default is 1}
\item{...}{Further arguments passed to \code{\link{update}}}
}

\value{
\item{rowpcoord}{Row principal coordinates}     
\item{colpcoord}{Column principal coordinates} 
\item{rowcoord}{Row standard coordinates}     
\item{colcoord}{Column standard coordinates} 
\item{sv}{Singular values} 
\item{inertia.e}{Percentages of explained (adjusted) inertia}
\item{levelnames}{Attribute names} 
\item{rowctr}{Row contributions}
\item{colctr}{Column contributions} 
\item{rowcor}{Row squared correlations} 
\item{colcor}{Column squared correlations} 
\item{rowmass}{Row masses} 
\item{colmass}{Column masses} 
\item{indmat}{Indicator matrix} 
\item{m}{Number of cases processed up to this point} 
\item{ff}{A copy of \code{ff} in the return object}
}

\references{
Ross, D. A., Lim, J., Lin, R. S., & Yang, M. H. (2008). Incremental learning for robust visual tracking. \emph{International Journal of Computer Vision}, \emph{77}(1-3), 125-141.\cr \cr
Iodice D' Enza, A., & Markos, A. (2015). Low-dimensional tracking of association structures in categorical data, \emph{Statistics and Computing}, \emph{25}(5), 1009-1022. 
}

\seealso{
 \code{\link{add_es}}, \code{\link{i_mca}}, \code{\link{plot.i_mca}}
}

\examples{
data(women, package = "idm")
dat = women[,c(1:4)]
res_mca = i_mca(dat[1:300,])
nchunk = seq(301,2107,258)
for (k in c(1:(length(nchunk)-1)))
{
  res_mca = update(res_mca,dat[c((nchunk[k]+1):nchunk[k+1]),],nchunk=k+1)
}
plot(res_mca, what=c(FALSE, TRUE), animation = FALSE)
}
