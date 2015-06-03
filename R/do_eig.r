do_eig <- function(data) {
 # require("corpcor")
  # data: data matrix
  out = list()
  n = nrow(data)  ## number of rows
  p = ncol(data)  ## number of columns
  
  #if (is.cate == F) {
  orgn = apply(data, 2, mean)
  orgn = t(t(orgn))
  # if (is.svd == T) {
  data = t(data)
  #  }
  oner = matrix(1, 1, n)
  #     if (is.svd == F) {
  #       cov.data = (1/n) * (t(data - t(oner) %*% t(orgn)) %*% (data - t(oner) %*% t(orgn)))
  #       eig.res = eigen(cov.data)
  #       #   print('cov.data')
  #       #   print(cov.data)
  #       vct = eig.res$vectors
  #       val = eig.res$values
  #       val = t(t(val))
  #     #  out = list()
  #       out$cov = cov.data
  #     } else {
  #       
  cen.mat = data - (orgn) %*% (oner)
  svd.res = fast.svd(cen.mat,0)
  vct = svd.res$u
  val = svd.res$d
  vctCol = svd.res$v
  #  out = list()
  out$vctCol = vctCol
  #out$cen.mat = cen.mat
  #    }
  #out=list()
  out$n = n
  out$orgn = orgn
  out$vct = vct
  out$val = val
  out
  #  }
  
  
  #   if (is.cate == T) {
  #     #######################
  # #     disdata = disjMake(data)  # disjunctive coding
  # #     J = disdata$J  #  number of modalities
  # #     Q = disdata$Q  #  number of variables
  # #     dZ = disdata$dZ  ## indicator matrix
  # dZ=as.matrix(dummy.data.frame(data,drop=F))
  # J=ncol(dZ)
  # Q=ncol(data)
  # #######################
  #     
  #     c = apply(dZ/(Q * (n)), 2, sum)  # column margins
  #     D = as.vector(c)
  #     SqD = as.vector(sqrt(D))
  #     invSqD= as.vector(1/sqrt(D))
  #     oner = matrix(1, n, 1)
  #     onec = matrix(1, J, 1)
  #     orgn = t(t((1/sqrt(n)) %*% t(onec)) * SqD)
  #     orgn = (t(orgn))  # mean vector
  #     Zp = t(t((dZ/(Q * sqrt(n)))) * invSqD)
  #     Zc = Zp - oner %*% t(orgn)  # cented version of the matrix
  #     
  #     if (is.svd == F) {
  #       
  #       Sz = t(Zc) %*% Zc  ## standardized residual matrix
  #       eig.Sz = eigen(Sz)  #EVD
  #       #eig.Sz = eigen(Sz*4)  #EVD
  #       # print(eig.Sz)
  #       vct = eig.Sz$vectors  # eigenvectors
  #       val = eig.Sz$values  # eigenvalues
  #       val = t(t(val))
  #       out$Sz = Sz
  #     }  #
  #     if (is.svd == T) {
  #       #print('inside')
  #       #svd.Zc = fast.svd(Zc)
  #       svd.Zc = svd(Zc)
  #       vct = svd.Zc$v
  #       vctCol = svd.Zc$u
  #       val = svd.Zc$d
  #       out$vctCol = vctCol
  #     }
  #     out$dZ = dZ
  #     out$Zc = Zc
  #     out$Zp = Zp
  #     out$cm = c
  #     out$Q = Q
  #     out$invSqD = invSqD
  #     out$D = D
  #     out$orgn = orgn
  #     out$n = n
  #     #out$Sz = Sz
  #     out$vct = vct
  #     out$val = val
  #     
  #     out
  #     
  #   }
  #   
  out
  
} 