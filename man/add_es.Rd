\name{add_es}
\alias{add_es}
\title{
Adds two eigenspaces using block-wise incremental SVD (with or without mean update)
} 
\description{
This function implements two procedures for updating existing decomposition. When \code{method="esm"} it adds two eigenspaces using the incremental method of Hall, Marshall & Martin (2002). The results correspond to the eigenspace of the mean-centered and concatenated data.
When \code{method = "isvd"} it adds the eigenspace of an incoming data block to an existing eigenspace using the block-wise incremental singular value decomposition (SVD) method described by Zha & Simon (1999), Levy and Lindenbaum (2000), Brand (2002) and Baker (2012). New data blocks are added row-wise. The procedure can optionally keep track of the data mean using the orgn argument, as described in Ross et al. (2008) and Iodice D'Enza & Markos (2015).
}
\usage{
add_es(eg, eg2, current_rank, ff = 0, method = c("esm", "isvd"))
}

\arguments{
\item{eg}{A list describing the eigenspace of a data matrix, with components \cr
\code{u} {Left eigenvectors} \cr
\code{v} {Right eigenvectors} \cr
\code{m} {Number of cases} \cr
\code{d} {Eigenvalues} \cr
\code{orgn} {Data mean} \cr
}
\item{method}{refers to the procedure being implemented: \code{"esm"} refers to the eigenspace merge (Hall et al., 2002); \code{"isvd"} refers to the incremental SVD method, with or without keeping track of the data mean.}
\item{eg2}{(*)A list describing the eigenspace of a data matrix, with components \cr
\code{u} {Left eigenvectors} \cr
\code{v} {Right eigenvectors} \cr
\code{m} {Number of cases} \cr
\code{d} {Eigenvalues} \cr
\code{orgn} {Data mean} \cr
}
\item{current_rank}{Rank of approximation; if empty, the full rank is used}
\item{ff}{(**)Number between 0 and 1 indicating the forgetting factor used to down-weight the contribution of earlier data blocks to the current solution. When ff = 0 (default) no forgetting occurs}
(*) for \code{method = "esm"} only; (**) for \code{method = "isvd"} only.
}

\value{
A list describing the SVD of a data matrix, with components \cr
\item{u}{Left singular vectors}
\item{d}{Singular values}
\item{v}{Right singular vectors}
\item{m}{Number of cases}
\item{orgn}{Data mean; returned only if \code{orgn} is given as input}
}
\references{
Zha, H., & Simon, H. D. (1999). On updating problems in latent semantic indexing. \emph{SIAM Journal on Scientific Computing}, \emph{21}(2), 782-791.\cr \cr
Levy, A., & Lindenbaum, M. (2000). Sequential Karhunen-Loeve basis extraction and its application to images. \emph{IEEE Transactions on Image Processing}, \emph{9}(8), 1371-1374.\cr \cr
Brand, M. (2002). Incremental singular value decomposition of uncertain data with missing values. In \emph{Computer Vision-ECCV 2002} (pp. 707-720). Springer Berlin Heidelberg.\cr \cr
Ross, D. A., Lim, J., Lin, R. S., & Yang, M. H. (2008). Incremental learning for robust visual tracking. \emph{International Journal of Computer Vision}, \emph{77}(1-3), 125-141.\cr \cr
Baker, C. G., Gallivan, K. A., & Van Dooren, P. (2012). Low-rank incremental methods for computing dominant singular subspaces. \emph{Linear Algebra and its Applications}, \emph{436}(8), 2866-2888.\cr \cr
Iodice D' Enza, A., & Markos, A. (2015). Low-dimensional tracking of association structures in categorical data, \emph{Statistics and Computing}, \emph{25}(5), 1009-1022. 
}

\seealso{
\code{\link{do_es}}, \code{\link{i_pca}}, \code{\link{i_mca}}, \code{\link{update.i_pca}}, \code{\link{update.i_mca}}}

\examples{

## Example 1 - eigenspace merge (Hall et al., 2002)
#Iris species
data("iris", package = "datasets")
X = iris[,-5]
#obtain two eigenspaces
eg = do_es(X[1:50, ])
eg2 = do_es(X[c(51:150), ])
#add the two eigenspaces keeping track of the data mean
eg12 = add_es(method = "esm", eg, eg2)
#equivalent to the SVD of the mean-centered data (svd(scale(X, center = TRUE,scale = FALSE)))

## Example 2 - block-wise incremental SVD with mean update, full rank (Ross et al., 2008)
data("iris", package = "datasets")
# obtain the eigenspace of the first 50 Iris species
X = iris[,-5]
eg = do_es(X[1:50, ])
#update the eigenspace of the remaining species to
eg_new = add_es(method = "isvd", eg, data.matrix(X[c(51:150), ]))
#equivalent to the SVD of the mean-centered data (svd(scale(X, center = TRUE, scale = FALSE)))

##Example 3 - incremental SVD with mean update, 2d approximation (Ross et al., 2008)
data("iris", package = "datasets")
# obtain the eigenspace of the first 50 Iris species
X = iris[,-5]
eg = do_es(X[1:50, ])
#update the eigenspace of the remaining species to
eg = add_es(method = "isvd", eg, data.matrix(X[c(51:150), ]),current_rank = 2)
#similar to PCA on the covariance matrix of X (SVD of the mean-centered data)

}

