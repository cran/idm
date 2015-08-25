\name{add_svd}
\alias{add_svd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adds two eigenspaces using incremental SVD (with or without mean update)
} 
\description{
This function adds the eigenspace of an incoming data block to an existing eigenspace using the incremental SVD method described by Zha & Simon (1999), Levy and Lindenbaum (2000), Brand (2002) and Baker (2012). New data blocks are added row-wise. The procedure can optionally keep track of the data mean using the orgn argument, as described in Ross et al. (2008) and Iodice D'Enza & Markos (2015).
}
\usage{
add_svd(eg,B,m,current_rank,orgn,ff = 0)
}

\arguments{
\item{eg}{A list describing the SVD of a data matrix, with components \cr
\code{u} {Left singular vectors} \cr
\code{v} {Right singular vectors} \cr
\code{d} {Singular values} \cr
}
\item{B}{Incoming data block (rows > columns)}
\item{m}{Number of cases (rows) processed up to this point}
\item{current_rank}{Rank of approximation; if empty then full rank is used}
\item{orgn}{Data mean; if empty then data is assumed as mean-centered}
\item{ff}{Number between 0 and 1 indicating the forgetting factor used to down-weight the contribution of earlier data blocks to the current solution. When ff = 0 (default) no forgetting occurs}
}

\value{
A list describing the SVD of a data matrix, with components \cr
\item{u}{Left singular vectors}
\item{d}{Singular values}
\item{v}{Right singular vectors}
\item{m}{Number of cases}
\item{orgn}{Data mean; returned only if 'orgn' is given as input}
}
\references{
Zha, H., & Simon, H. D. (1999). On updating problems in latent semantic indexing. \emph{SIAM Journal on Scientific Computing}, 21(2), 782-791.\cr \cr
Levy, A., & Lindenbaum, M. (2000). Sequential Karhunen-Loeve basis extraction and its application to images. \emph{IEEE Transactions on Image Processing}, 9(8), 1371-1374.\cr \cr
Brand, M. (2002). Incremental singular value decomposition of uncertain data with missing values. In \emph{Computer Vision-ECCV 2002} (pp. 707-720). Springer Berlin Heidelberg.\cr \cr
Ross, D. A., Lim, J., Lin, R. S., & Yang, M. H. (2008). Incremental learning for robust visual tracking. \emph{International Journal of Computer Vision}, 77(1-3), 125-141.\cr \cr
Baker, C. G., Gallivan, K. A., & Van Dooren, P. (2012). Low-rank incremental methods for computing dominant singular subspaces. \emph{Linear Algebra and its Applications}, 436(8), 2866-2888.\cr \cr
Iodice D' Enza, A., & Markos, A. (2015). Low-dimensional tracking of association structures in categorical data, \emph{Statistics and Computing}, 25(5), 1009-1022. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{add_eig}}, \code{\link{update_mca}}, \code{\link{i_pca}}, \code{\link{i_mca}}
}
\examples{
##Example 1 - Update SVD (full rank)
data("iris", package = "datasets")
#obtain the SVD of the first 50 Iris species
X = iris[,-5]
eg = svd(X[1:50,])
#update the eigenspace of the remaining species to
eg = add_svd(eg,data.matrix(X[c(51:150),]),50)
#similar results to svd(X)

##Example 2 - Update SVD with mean update (full rank)
data("iris", package = "datasets")
#obtain the eigenspace of the first 50 Iris species
X = iris[,-5]
eg = do_eig(X[1:50,])
#update the eigenspace of the remaining species to
eg = add_svd(eg,data.matrix(X[c(51:150),]),50,,eg$orgn)
#eg$v corresponds to PCA loadings on the covariance matrix of X

##Example 3 - Update SVD (low-rank approximation: current_rank = 2)
data("iris", package = "datasets")
#obtain the SVD of the first 50 Iris species
X = iris[,-5]
eg = svd(X[1:50,])
#update the eigenspace of the remaining species to
eg = add_svd(eg,data.matrix(X[c(51:150),]),50,2)
#similar results to svd(X)
}

