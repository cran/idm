update.i_pca <- function(object, incdata, current_rank, ...) {
  #This is equivalent to the PCA on the covariance (checked 04 Feb 2015)
  #Fixed and updated 23 Mar 2017
  #data <- data.frame(lapply(data.frame(incdata), factor))
  
  eg = object$eg
  data = data.matrix(incdata)
  
  if (length(object$levelnames) != dim(incdata)[2]){stop("The data blocks must have the same number of columns/variables")}
 
  if (missing("current_rank")) {
    #current_rank = full rank
    current_rank = dim(incdata)[2]
  }
  
  nrows = eg$m + nrow(data)
  ncols = ncol(data)
  collabs = colnames(incdata)
  dims = current_rank #ncols
  #keep these just for sign check
  PC1 = eg$v[,1:current_rank]
  PCu1 = eg$u[,1:current_rank]%*%diag(eg$d[1:current_rank]) 
 
  mat.chu = data
  #eigenspace of the incoming block
  eg2 = do_es(mat.chu)
  #add eigenspaces
  eg = add_es(eg, eg2,current_rank = current_rank, method="esm")
  PCall = (1/sqrt(nrows))*eg$v[,1:current_rank]%*%diag(eg$d[1:current_rank])
  PCuall = eg$u[,1:current_rank]%*%diag(eg$d[1:current_rank]) 
  nrows2 = nrow(mat.chu)    
  ## insert ctr comps
  signe = 2*(PCuall>0)-1
  PCuall2 = PCuall^2
  MF2 = (1/nrows2)*PCuall2
  V = apply(MF2,2,sum)
  dist2_12.var = apply(PCall^2,1,sum)
  # Contributions of observations to the components
  PCuall.ctr = MF2%*%diag(rep(1,current_rank)/V)#*signe
  PCall.ctr = t(t(PCall^2)/(eg$d * (1/sqrt(eg$m)))^2)*100
  # Squared cosines of the observations
  PCuall.cor = suppressWarnings(MF2 / (apply(MF2,1,sum)*rep(1,ncols)))
  #TODO: COLCOR (only) are not similar to i_pca and current_rank less than full and close to 2!
  PCall.cor = (PCall / sqrt(dist2_12.var))^2
  
  PCall = sign_match(PC1, PCall)
  PCuall = sign_match(PCu1, PCuall)
  out = list()
  out$eg=eg
#  out$u = eg$u
#  out$v =  eg$v
  
  # PCA eigenvalues
  sv = eg$d/sqrt(nrows)
  if (current_rank == ncols) {
    out$inertia.e=eg$d^2/(sum(eg$d^2)) #sv/(sum(sv)^2)
  } else {
    out$inertia.e= sv^2/ncols 
  }
  
  
  
#  out$d = eg$d 
  #  out$eg=eg
  out$m = nrows
  out$rowctr=PCuall.ctr[,c(1:dims)]*100
  out$rowcor = PCuall.cor[,c(1:dims)]
  out$colctr= PCall.ctr[,c(1:dims)]
  out$colcor = PCall.cor[,c(1:dims)]
  out$levelnames = collabs
  out$colpcoord = PCall
  out$sv = sv
  out$rowpcoord = PCuall 
  class(out)="i_pca"
  return(out)
}
