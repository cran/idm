update.i_mca <-  function(object, incdata, current_rank, ff = 0,  ...) {
  #Fixed and updated 23 Mar 2017
  eg = list()
  # ff = 1 - ff
  data <- data.frame(lapply(data.frame(incdata), factor))
  
  mods1 = apply(data, 2, unique)
  if(is.null(dim(mods1))== FALSE){
    mods1=split(t(mods1),colnames(mods1))
  } 
  
  nchunk=1 #check i_mca calculation of row/col principal coords to see why it's needed
  
  Q = ncol(incdata)
  n1 = object$m 
  if (missing("current_rank")) {
    #full rank
    current_rank =  length(object$sv)
  }
  
  labs = object$levelnames
  r = object$rowmass
  c = object$colmass
  n.mods1 = object$levels.n
  J = ncol(object$indmat)
  eg$orgn = object$orgn
  SRall = object$rowcoord[,1:current_rank] 
  SCall = object$colcoord[,1:current_rank] 
  eg$m = n1
  eg$u = SRall * sqrt(r)
  eg$v = SCall * sqrt(c)
  eg$d = object$sv
  
  dims = current_rank
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
  #center sZ2!!
  # sZ2 = sZ2 - rep(eg$orgn, rep.int(nrow(sZ2), ncol(sZ2)))
  ###### END OF THE NEW CODE ########################
  n2 = n.chu
  n12= n1 + n2
  c12 = (c*n1 + c2*n2)/n12
  r12 = rep(1/n12,n12)
  eg = add_es(eg,sZ2,current_rank,ff=ff,method="isvd")
  #update column standard coordinates
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
  PCall = sign_match(SRall, PCall[,1:current_rank])
  PCuall = sign_match(SCall, PCuall[,1:current_rank])
  
  PCall.ctr=eg$v^2
  PCall.cor=(PCall^2)/apply(PCall^2,1,sum)
  PCuall.ctr =  suppressWarnings(t(((1/n1)*t(PCuall^2))*as.vector(1/(eg$v)^2)))
  PCuall.cor = (PCuall^2)/apply(PCuall^2,1,sum) 
  
  #### Inertia Percentages, COR, CTR are approximate #05.02.2016
  
  #calculates  inertia
  nd.max = ncol(tZ2$SZ)-Q
  eg$d = eg$d/sqrt(nchunk+1)
  if (current_rank == nd.max) {
    inertia0    <- eg$d[1:nd.max]^2
    inertia.t   <- nd.max/Q
    inertia.e   <- inertia0 / inertia.t
  } else {
    inertia0    <- eg$d[1:current_rank]^2
    inertia.t   <- nd.max/Q
    inertia.e   <- inertia0 / inertia.t
  }
  
  out = list()
  out$rowpcoord = PCuall[,c(1:dims)]/sqrt(nchunk+1) #/sqrt(nchunk)
  out$colpcoord = PCall[,c(1:dims)]/sqrt(nchunk+1) #/sqrt(nchunk)
  out$rowcoord = SRall[,c(1:dims)]  
  out$colcoord = SCall[,c(1:dims)]
  out$sv = eg$d
  out$inertia.e = inertia.e #inertia.adj / inertia.t
  out$levelnames = labs
  out$rowctr = PCuall.ctr[,c(1:dims)]
  out$colctr = PCall.ctr[,c(1:dims)]
  out$rowcor = PCuall.cor[,c(1:dims)]
  out$colcor = PCall.cor[,c(1:dims)]
  out$rowmass = r
  out$colmass = c
  out$indmat = tZ2$dZ
  out$orgn = eg$orgn
  out$m = n1 
  out$ff = ff
  class(out)="i_mca"
  return(out)
}   
