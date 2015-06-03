update_mca <-  function(mca, incdata, f = 0) {
  eg = list()
  f = 1 - f
  data <- data.frame(lapply(data.frame(incdata), factor))
  
  mods1 = apply(data, 2, unique)
  if(is.null(dim(mods1))== FALSE){
    mods1=split(t(mods1),colnames(mods1))
  }    
  
  
  Q = ncol(incdata)
  #m = mca$m
  n1 = mca$m #nrow(Z)
  labs = mca$levelnames
  r = mca$rowmass
  c = mca$colmass
  n.mods1 = mca$levels.n
  J = ncol(mca$indmat)
  eg$orgn = mca$orgn
  SRall = mca$rowcoord
  SCall = mca$colcoord
  
  eg$m = n1
  eg$u = SRall * sqrt(r)
  eg$v = SCall * sqrt(c)
  eg$d = mca$sv
  dims=dim(eg$v)[2]
  current_rank = dim(eg$v)[2]  ## This is J-Q
  mat.chu = data 
  n.chu = nrow(mat.chu)
  
  mods.up = apply(mat.chu, 2, unique)
  if(is.null(dim(mods.up))== F){
    mods.up=split(t(mods.up),colnames(mods.up))
  }    
  
  n.mods.up = sapply(mods.up, length)
  catDiff = n.mods1 - n.mods.up
  
  if (any(catDiff) != 0) {
    fake.row = fake_row_make(mods1, mods.up, n.mods1, n.mods.up)
    mat.chu = rbind(mat.chu, fake.row)
    tZ2 = transform_z(mat.chu, is.weight = T, is.exact = F, c = c)#, r = r, c = c)
    sZ2 = tZ2$SZ[1:n.chu, ]
    c2 = tZ2$c
  } else {
    tZ2 = transform_z(mat.chu, is.weight = T, is.exact = F, c = c)#, r = r, c = c)
    sZ2 = tZ2$SZ[1:n.chu, ]
    c2 = tZ2$c
  }
  ###### END OF THE NEW CODE ########################
  n2 = n.chu
  n12= n1 + n2
  c12 = (c*n1 + c2*n2)/n12
  r12 = rep(1/n12,n12)
  eg = add_svd(eg$u, eg$d, eg$v, sZ2, eg$m, eg$orgn,f, current_rank)
  #  eg2 = do_eig(sZ2)
  #  eg$vctCol = eg$u
  #  eg$vct = eg$v
  #  eg$val = eg$d
  #  eg$n = eg$m
  #  eg = add_eig(eg,eg2)
  #  eg$u = eg$vctCol 
  #  eg$v =eg$vct 
  #  eg$d = eg$val 
  #  eg$m = eg$n 
  #update column standard coordinates
  # eg$d = eg$d/sqrt(2)
  #eg$v = sqrt(eg$v^2/2)
  SCall <- (eg$v / sqrt(c12))
  #update column principal coordinates
  PCall <- SCall %*% diag(as.vector((eg$d[1:current_rank])))
  #update row standard coordinates
  SRall <- (eg$u / sqrt(r12))
  #update row principal coordinates
  PCuall <- SRall[,c(1:current_rank)] %*% diag(as.vector((eg$d)))
  #PCuall <- Z %*% as.matrix(SCall) * Q^(-1)
  n1 = n12
  c = c12
  r = r12
  
  PCall = sign_match(mca$colpcoord, PCall[,1:current_rank])
  PCuall = sign_match(mca$rowpcoord, PCuall[,1:current_rank])
  
  PCall.ctr=eg$v^2
  PCall.cor=(PCall^2)/apply(PCall^2,1,sum)
  PCuall.ctr =  suppressWarnings(t(((1/n1)*t(PCuall^2))*as.vector(1/(eg$v)^2)))
  PCuall.cor = (PCuall^2)/apply(PCuall^2,1,sum) 
  
  #calculates adjusted inertia
  J = length(eg$d)
  #get (almost) same eigenvalues with exact
  dd = eg$d # sqrt(eg$d^2/2)
  inertia0 = dd[1:(J-Q)]^4
  alldim <- sum(sqrt(inertia0) >= 1/Q)
  inertia.adj  <- ((Q/(Q-1))^2 * (sqrt(inertia0)[1:alldim] - 1/Q)^2)
  inertia.t    <- (Q/(Q-1)) * (sum(inertia0) - ((J - Q) / Q^2))
  
  out = list()
  out$sv = eg$d
  #out$inertia_e = (eg$d[c(1:alldim)])^2/inertia.t
  out$inertia_e = inertia.adj / inertia.t
  out$rowcoord = SRall[,c(1:dims)]
  out$rowpcoord = PCuall[,c(1:dims)]
  # out$colpcoord = sqrt(SCall[,c(1:dims)]^2/2)
  out$colcoord = SCall[,c(1:dims)]
  out$colpcoord = PCall[,c(1:dims)]
  out$colctr = PCall.ctr[,c(1:dims)]
  out$rowctr = PCuall.ctr[,c(1:dims)]
  out$colcor = PCall.cor[,c(1:dims)]
  out$rowcor = PCuall.cor[,c(1:dims)]
  out$rowmass = r
  out$colmass = c
  out$indmat = tZ2$dZ
  out$levelnames = labs
  out$orgn = eg$orgn
  out$m = n1 
  out$f = 1-f
  class(out)="i_mca"
  return(out)
}   
