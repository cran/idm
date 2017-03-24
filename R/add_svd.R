add_svd <- function(eg,B,m,current_rank,ff = 0) {
  out = list()
  #new data block
  B = data.matrix(B) 
  #columns (fixed)
  c = dim(B)[2]
  #num of new rows
  n = dim(B)[1]
  if (missing("current_rank")) {
    #full rank
    current_rank = c
  }
  #for convenience
  ff = 1 -ff
  #get low-rank matrices
  Pk = eg$u[,1:current_rank]
  Sk = eg$d[1:current_rank]
  Qk = eg$v[,1:current_rank]
  
  
  #mean calculation
  orgn = eg$orgn
  orgnb = colMeans(B) 
  #center data
  Bc = scale(B,center=TRUE,scale=FALSE)
  #account for the variance of the mean
  Bc <- rbind(Bc,t(sqrt((n*m)/(n+m))*as.matrix((orgnb-orgn))))
  #mean update
  orgnc <- (ff*m*orgn + n*orgnb)/(n+ff*m)
  out$orgn <- orgnc
  #QR-decomposition of (I-Qk'Qk)B
  
  qrstr = qr((diag(c) - Qk%*%t(Qk))%*%t(Bc))
  Qt = qr.Q(qrstr,complete=T)
  L = qr.R(qrstr,complete=T)
  
  #Form middle matrix K
  #fix this when current_rank=1
  K = rbind(cbind(ff*diag(Sk),matrix(0,current_rank,c)),cbind(Bc%*%Qk,t(L)))
  #get svd of K
  eg = fast.svd(K,0)
  #keep current_rank singular values and vectors
  Uk = eg$u[,1:current_rank]
  Sk = eg$d[1:current_rank]
  Vk = eg$v[,1:current_rank]
  #calculate U and V
  Vk = cbind(Qk,Qt)%*%Vk
  
  #Fixed 24.02.17 during the industrial affiliates conference
  # Uk = rbind(cbind(Pk,matrix(0,m,n+1)),cbind(matrix(0,n,current_rank+1),diag(n)))%*%Uk
  Uk = rbind(cbind(Pk,matrix(0,m,n+1)),cbind(matrix(0,n,current_rank),cbind(diag(n),rep(0,n,1))))%*%Uk
  #update number of columns m processed so far
  m <- m + n
  #but the correct is 
  #  m <- ff*m + n
  
  #output
  out$u <- Uk
  out$d <- Sk
  out$v <- Vk
  out$m <- m 
  
  
  out
  
}

