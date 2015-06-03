update_pca <- function(eg, incdata) {
  #eg = list()
  #data <- data.frame(lapply(data.frame(incdata), factor))
  data = incdata
  nrows = eg$n + nrow(data)
  ncols = ncol(data)
  collabs = colnames(incdata)
  dims = ncols
  PC1 = eg$vct
  PCu1 = eg$vctCol%*%diag(eg$val) 
  
  mat.chu = data
  eg2 = do_eig(mat.chu)
 
  eg12 = add_eig(eg, eg2)
  PCall = eg12$vct
  PCuall = eg12$vctCol%*%diag(eg12$val) 
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
  
  eg = eg12
  
  out = list()
  out$vctCol =  eg$vctCol #PCuall[,c(1:dims)]  
  out$vct =  PCall[,c(1:dims)]
  
  # PCA eigenvalues
  sv = eg12$val#/sqrt(nrows)
  # out$inertia_e=sv/(sum(sv))
  out$val = sv[c(1:dims)] 
  out$inertia_e=sv/(sum(sv))
  out$n = nrows
  out$rowctr=PCuall.ctr[,c(1:dims)]
  #rownames(out$rowctr) = rowlabs
  out$rowcor=PCuall.cor[,c(1:dims)]
  out$levelnames = collabs
  out$orgn = eg12$orgn
  out$colpcoords = out$vct
  out$rowpcoords = PCuall
  class(out)="i_pca"
  return(out)
}
