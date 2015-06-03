add_eig <- function(m1, m2) {
  
  #output
  out = list()
  
  bign = m1$n + m2$n  
  bigOrgn = ((m1$orgn * m1$n) + (m2$orgn * m2$n))/bign
  
  dorgn = m1$orgn - m2$orgn
  G = t(m1$vct) %*% m2$vct
  H = m2$vct - m1$vct %*% G
  g = t(m1$vct) %*% dorgn
  h = dorgn - m1$vct %*% g
  eps = 2.2204e-16
  HH = apply(H * H, 2, sum)
  hh = apply(h * h, 2, sum)
  if (any(HH > eps) | any(hh > eps)) {
 #   source("orth.r")
    Hnz = which(HH > eps)
    hnz = which(hh > eps)
    nu = orth(cbind(HH[, Hnz], hh[, hnz]))
    H = t(m1$vct) %*% nu
    nu = nu[, -Hnz]
  } else {
    nu = matrix(0, length(m1$orgn), 0)
  }
  resn1 = ncol(m1$vct) - nrow(m1$vct)
#  if (resn1 > 0) {
#    rpern1 = m1.resudue/resn1
#  } else {
    rpern1 = 0
#  }
  resn2 = ncol(m2$vct) - nrow(m2$vct)
#  if (resn2 > 0) {
#    rpern2 = m2.resudue/resn2
#  } else {
    rpern2 = 0
#  }
  nnu = nrow(nu)
  mnu = ncol(nu)
  Gamma = t(nu) %*% m2$vct
  
  D = t(t(G) * as.vector(m2$val))
  E = t(t(Gamma) * as.vector(m2$val))
  gamma = t(nu) %*% dorgn
#   if (is.svd == F) {
#     chk1 = cbind(t(m1$val), rpern1 %*% matrix(1, 1, mnu))
#     A1 = (m1$n/bign) * diag(as.vector(chk1))
#     #A1try = (m1$n/bign) * (as.vector(chk1))
#     A21 = cbind(D %*% t(G), D %*% t(Gamma))
#     A22 = cbind(E %*% t(G), E %*% t(Gamma))
#     A2 = (m2$n/bign) * rbind(A21, A22)
#     A31 = cbind(g %*% t(g), g %*% t(gamma))
#     A32 = cbind(gamma %*% t(g), gamma %*% t(gamma))
#     A3 = ((as.double(m1$n) * as.double(m2$n))/(bign^2)) * rbind(A31, A32)
#     A = A1 + A2 + A3
#     A = (A + t(A))/2
#     eigA = eigen(A)
#     vct = eigA$vectors
#     vct = cbind(m1$vct, nu) %*% vct
#     val = eigA$values
#   } else {
    A11 = cbind(m1$val * t(m1$vctCol), G %*% (m2$val * t(m2$vctCol)) )
    A12b = (t(nu) %*% m2$vct) %*% (m2$val * t(m2$vctCol))
    t = nrow(A12b)
    A12a = matrix(0, t, (m1$n))
    A12 = cbind(A12a, A12b)
    A1 = rbind(A11, A12)
    onerX = matrix(1, 1, m1$n)
    onerY = matrix(1, 1, m2$n)
    #dorgnX = m1$orgn - bigOrgn
    dorgnX = as.double(m1$orgn) - as.double(bigOrgn)
    dorgnY = as.double(m2$orgn) - as.double(bigOrgn)
    A21a = (t(m1$vct) %*% as.double(dorgnX) %*% onerX)
    A21b = (t(m1$vct) %*% as.double(dorgnY) %*% onerY)
    #A21b=(t(m2$vct) %*% dorgnX %*% onerX) prova con Uy
    A21 = cbind(A21a, A21b)
    A22a = t(nu) %*% as.double(dorgnX) %*% onerX
    A22b = t(nu) %*% as.double(dorgnY) %*% onerY
    A22 = cbind(A22a, A22b)
    A2 = rbind(A21, A22)
    A = A1 + A2
    #library(svd)
    svd.res = fast.svd(A,0)
    
 #   svd.res = irlba(A,nrow(A),nrow(A))
    #svd.res = svd(A)
    vct = svd.res$u
    vct = cbind(m1$vct, nu) %*% vct
    #vct=sign.matcher(m1$vct,vct)
    val = svd.res$d
    vctCol = svd.res$v
    out$vctCol = vctCol
#  }
  
  out$n = bign
  out$orgn = bigOrgn
  out$vct = vct
  out$val = val
 # out$A = A
  out
}  

