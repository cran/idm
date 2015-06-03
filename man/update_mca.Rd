\name{update_mca}
\alias{update_mca}
\title{
Updates a Multiple Correspondence Analysis solution
}
\description{
This function updates the Multiple Correspondence Analysis (MCA) solution on the indicator matrix using the incremental method of Ross et al. (2008)
}
\usage{
update_mca(mca, incdata, f = 0)
}

\arguments{
  \item{incdata}{
  Matrix of incoming data
}
\item{mca}{
A list describing an MCA solution with components \cr
\code{indmat} {Indicator matrix} \cr
\code{colcoord} {Column standard coordinates} \cr
\code{rowcoord} {Row standard coordinates} \cr    
\code{colpcoord} {Column principal coordinates} \cr
\code{rowpcoord} {Row principal coordinates}\cr     
\code{colctr} {Column contributions} \cr
\code{colcor} {Column squared correlations} \cr
\code{rowctr} {Row contributions} \cr
\code{rowcor} {Row squared correlations} \cr
\code{sv} {Eigenvalues} \cr
\code{rowmass} {Row masses} \cr
\code{colmass} {Column masses} \cr 
\code{levelnames} {Attribute names} \cr
\code{m} {Number of cases processed up to this point} \cr
\code{orgn} {Data mean} \cr
\code{inertia.e} {Percentages of explained (adjusted) inertia}}
\item{f}{Number between 0 and 1 indicating the "forgetting factor" used to down-weight the contribution of earlier data blocks to the current solution. When f = 0 (default) no forgetting occurs}
}

\value{
A list describing an MCA solution with components 
\item{indmat}{Indicator matrix} 
\item{colcoord}{Column standard coordinates} 
\item{rowcoord}{Row standard coordinates}     
\item{colpcoord}{Column principal coordinates} 
\item{rowpcoord}{Row principal coordinates}     
\item{colctr}{Column contributions} 
\item{colcor}{Column squared correlations} 
\item{rowctr}{Row contributions}
\item{rowcor}{Row squared correlations} 
\item{sv}{Eigenvalues} 
\item{rowmass}{Row masses} 
\item{colmass}{Column masses} 
\item{levelnames}{Attribute names} 
\item{m}{Number of cases processed up to this point} 
\item{orgn}{Data mean} 
\item{inertia_e}{Percentages of explained (adjusted) inertia} 
\item{f}{Number between 0 and 1 indicating the "forgetting factor"}
}

\references{
Ross, D. A., Lim, J., Lin, R. S., & Yang, M. H. (2008). Incremental learning for robust visual tracking. \emph{International Journal of Computer Vision}, 77(1-3), 125-141.\cr
Iodice D' Enza, A., & Markos, A. (2014). Low-dimensional tracking of association structures in categorical data. \emph{Statistics and Computing}, (forthcoming).
}

\seealso{
\code{\link{update_pca}}, \code{\link{i_pca}}, \code{\link{add_svd}}, \code{\link{add_eig}}
}

\examples{
data(women, package="idm")
dat = women[,c(1:4)]
n0 = 301
outmca = mjca(dat[1:n0,],lambda="indicator",ret=TRUE)
outZ = transform_z(dat[1:n0,],is.weight=FALSE)
outmca$m = n0
outmca$rowmass = outZ$r
outmca$orgn = colMeans(outZ$SZ[1:n0,])
nchunk = seq(n0,2107,258)
for (k in c(1:(length(nchunk)-1)))
{
  outmca = update_mca(outmca,dat[c((nchunk[k]+1):nchunk[k+1]),])
}
#to get coordinates on a similar scale to the ones produced with batch mca
outmca$colpcoord = outmca$colpcoord/sqrt(length(nchunk))
outmca$rowpcoord = outmca$rowpcoord/sqrt(length(nchunk))

plot(outmca, what=c(FALSE, TRUE), animation = FALSE)
}
