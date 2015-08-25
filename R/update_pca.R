update_pca <- function(eg, incdata) {
  #This is equivalent to the PCA on the covariance (checked 23 Aug 2015)
  #data <- data.frame(lapply(data.frame(incdata), factor))
  data = data.matrix(incdata)
  nrows = eg$m + nrow(data)
  ncols = ncol(data)
  collabs = colnames(incdata)
  dims = ncols
  #keep these just for sign check
  PC1 = eg$v
  PCu1 = eg$u%*%diag(eg$d) 
  
  mat.chu = data
  #eigenspace of the incoming block
  eg2 = do_eig(mat.chu)
  #add eigenspaces
  eg = add_eig(eg, eg2)
  PCall = eg$v
  PCuall = eg$u%*%diag(eg$d) 
  nrows2 = nrow(mat.chu)    
  
  ## insert ctr comps
  signe = 2*(PCuall>0)-1
  PCuall2 = PCuall^2
  MF2 = (1/nrows2)*PCuall2
  V = apply(MF2,2,sum)
  # Contributions of observations to the components
  PCuall.ctr = MF2%*%diag(rep(1,ncols)/V)#*signe
  # Squared distance to the origin
  d2 = apply(MF2,1,sum)
  # Squared cosines of the observations
  PCuall.cor = suppressWarnings(MF2 / (d2*rep(1,ncols)))
  
  PCall = sign_match(PC1, PCall)
  PCuall = sign_match(PCu1, PCuall)
  
 # eg = eg12
  
  out = list()
  out$u = eg$u#  PCuall[,c(1:dims)]  #eg$u #
  out$v =  eg$v #PCall[,c(1:dims)]
  
  # PCA eigenvalues
  sv = eg$d/sqrt(nrows)
  out$inertia_e=sv/(sum(sv))
  out$d = eg$d 
  out$m = nrows
  out$rowctr=PCuall.ctr[,c(1:dims)]
  out$rowcor=PCuall.cor[,c(1:dims)]
  out$levelnames = collabs
  out$orgn = eg$orgn
  out$colpcoords = out$v
  out$sv = sv
  out$rowpcoords = PCuall 
  class(out)="i_pca"
  return(out)
}
