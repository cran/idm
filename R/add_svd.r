add_svd <- function(U,S,V,B,m,orgn,f,current_rank) {
  # add_svd <- function(U,S,V,B,m,orgn,f,current_rank) {
  #   Given the SVD of
  #   A = U*S*V' #and the orgn of A
  #   update it to be the SVD of
  #   [A ; B] = Up*Sp*Vp'
  #   that is, add new rows
  #   Currently, Up is accurate but Vp is approximated 
  #   This is equivalent to Brand's approach + mean update 
  #   which is described in Ross et al. (2007)
  
  # The subspace rotations involved may not preserve orthogonality due
  # to numerical round-off errors.  To compensate, you can set the
  # "force_orth" flag, which will force orthogonality via a QR plus
  # another SVD.  In a long loop, you may want to force orthogonality
  # every so often.
  out=list()
  # center B
  n = dim(B)[1]
  orgnb = colMeans(B) 
  Bc <- B - t(as.matrix(orgnb) %*% as.matrix(t(rep(1,n))))
  #account for the variance of the mean
  Bc <- rbind(Bc,t(sqrt((n*m)/(n+m))*as.matrix((orgnb-orgn))))
  #  orgnc <- (m / (m + n))*orgn + (n / (m + n))*orgnb
  orgnc <- (f*m / (f*m + n))*orgn + (n / (f*m + n))*orgnb
  # P is an orthogonal basis of the column-space
  # of (I-UU')a, which is the component of "a" that is
  # orthogonal to U.
  L <-  Bc %*% V
  p <- Bc - L %*% t(V)
  P <- orth(p)
  # p may not have full rank.  If not, P will be too small.  Pad
  # with zeros.
  P <- cbind(P,matrix(0,dim(P)[1],dim(p)[2]-dim(P)[2]))
  
  Ra <-  p %*% t(P)
  # Diagonalize K, maintaining rank
  b <- dim(diag(S))[2]
  #account for the variance of the mean  
  K <- cbind(f*rbind(diag(S),L),rbind(matrix(0, b, n+1),Ra))
  # K <- cbind(rbind(diag(S),L),rbind(matrix(0, b, n),Ra))
#  egk <- irlba(K,current_rank,current_rank)
  egk <- fast.svd(K,0)
  egk$u <- egk$u[,1:current_rank]
  egk$d <- egk$d[1:current_rank]
  egk$v <- egk$v[,1:current_rank]
  # Now update our matrices!
  Sp <- egk$d
  #account for the variance of the mean
  AA = cbind(rbind(U,matrix(0,n,current_rank)),rbind(matrix(0,m,n+1),cbind(diag(n),rep(0,n))))
  #AA = cbind(rbind(U,matrix(0,n,current_rank)),rbind(matrix(0,m,n),diag(n)))
  Up <- AA%*%egk$u
  # Up <- rbind(U, L)
  # Exploit structure to compute this fast: Vp = [ V Q ] * tVp;
  Vp = t(egk$v)%*%rbind(t(V),P)
  Vp <- t(Vp)[,1:current_rank]
  Sp <- Sp[1:current_rank]
  #Up <- rbind(U,P[,c(1:current_rank)]) %*% t(egk$u) 
  Up <- Up[,1:current_rank]
  # The above rotations may not preserve orthogonality, so we explicitly
  # deal with that via a QR plus another SVD.  In a long loop, you may
  # want to force orthogonality every so often.
  
  #Force orthogonality
 # URQ= qr(Up)
#  UQ = qr.Q(URQ)
#  UR = qr.R(URQ)
#  VRQ= qr(Vp)
#  VQ = qr.Q(VRQ)
#  VR = qr.R(VRQ)
#  egk <- fast.svd(UR %*% diag(Sp) %*% t(VR))
#  Up = UQ %*% egk$u
#  Vp = VQ %*% egk$v
#  Sp = egk$d
  #update m
  m <- m + n
  
  out$u <- Up
  out$d <- Sp
  out$v <- Vp
  out$m <- m 
  out$orgn <- orgnc
  out
  
}