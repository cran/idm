i_pca <- function(data1, data2, nchunk = 2, disk=FALSE) {
  ## data1 the starting matrix; data2 the upcoming data
  ## nchunk: number of chunks to split data2
  ## disk: store output to use for graphics
  
  ## This is equivalent to PCA on the covariance matrix 
  ## (i.e. on mean centered data)
  ##  princomp(data,cor=F) 
  # source("do_eig.r")
  #  source("add_eig.r")
  #  source("mat_split.r")
  #  source("sign_match.r")
  #  require("corpcor")
  
if(disk==TRUE){
  suppressWarnings(dir.create("./PCstory"))
}else{
  allCoords=list()
  allCoordsU=list()
  allctrU=list()
  allcorU=list()
}

rowlabs = c(rownames(data1),rownames(data2))
collabs = colnames(data1)
nrows = nrow(rbind(data1,data2))
ncols = ncol(data1)
nrows1 = nrow(data1)

dims=ncols

if ((length(nchunk) > 1 ) & (sum(nchunk) != nrow(data2))) {
  stop("\nchunk blocks do not match the size of 'data2'")
}


eg1 = do_eig(data1)
PC1 = eg1$vct
PCu1 = eg1$vctCol%*%diag(eg1$val) 

## insert ctr comps
signe = 2*(PCu1>0)-1
PCu12 = PCu1^2
MF2 = (1/nrows1)*PCu12
V = apply(MF2,2,sum)
# Contributions of observations to the components
PCu1.ctr = MF2%*%diag(rep(1,ncols)/V)#*signe
# Squared distance to the origin
d2 = apply(MF2,1,sum)
# Squared cosines of the observations
PCu1.cor = suppressWarnings(MF2 / (d2*rep(1,ncols)))

if(disk==TRUE){
  fnameA=paste("./PCstory/PCstart",1,".txt",sep="")
  fnameB=paste("./PCstory/PCEnd",1,".txt",sep="")
  fnameC=paste("./PCstory/PCstartUnit",1,".txt",sep="")
  fnameD=paste("./PCstory/PCendUnit",1,".txt",sep="")
  fnameE=paste("./PCstory/PCctrUnit",1,".txt",sep="")
  fnameF=paste("./PCstory/PCcorUnit",1,".txt",sep="")
  #   write.table(file=fnameA, matrix(0,dim(PC1[,1:dims])[1],dims))
  write.table(file=fnameA, PC1[,1:dims])
  write.table(file=fnameB, PC1[,1:dims])
  #   write.table(file=fnameC, matrix(0,dim(PCu1[,1:dims])[1],dims))
  write.table(file=fnameC,PCu1[,1:dims])   
  write.table(file=fnameD, PCu1[,1:dims])
  write.table(file=fnameE, PCu1.ctr[,1:dims])
  write.table(file=fnameF, PCu1.cor[,1:dims])
}

if(disk==FALSE){
  allCoordsU[[1]]=PCu1[,c(1:dims)]
  allCoords[[1]]=PC1[,c(1:dims)]
  allctrU[[1]] = PCu1.ctr[,c(1:dims)]
  allcorU[[1]] = PCu1.cor[,c(1:dims)]
}

out.split = mat_split(data2, (nchunk))  
mat.story = out.split$splitMat

#if block sizes are given, switch back to number
if (length(nchunk) > 1) {
  nchunk = length(nchunk)
} 

for (q in 1:length(mat.story)) {
  mat.chu = data.matrix(mat.story[[q]])
  
  ### coordinate computation
  ## column computation (modalities)
  if (q > 1) {
    #######################
    PCu1 = PCuall
    PC1 = PCall
    #######################
  }
  
  eg2 = do_eig(mat.chu)
  eg12 = add_eig(eg1, eg2)
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
  
  if(disk==FALSE){      
    allCoords[[q+1]]=PCall[,c(1:dims)]
    allCoordsU[[q+1]]=PCuall[,c(1:dims)]
    allctrU[[q+1]] = PCuall.ctr[,c(1:dims)]
    allcorU[[q+1]] = PCuall.cor[,c(1:dims)]
  } 
  
  eg1 = eg12
  
  if(disk==TRUE){      
    fnameA=paste("./PCstory/PCstart",q,".txt",sep="")
    fnameB=paste("./PCstory/PCEnd",q+1,".txt",sep="")
    fnameC=paste("./PCstory/PCstartUnit",q,".txt",sep="")
    fnameD=paste("./PCstory/PCendUnit",q+1,".txt",sep="")
    fnameE=paste("./PCstory/PCctrUnit",q+1,".txt",sep="")
    fnameF=paste("./PCstory/PCcorUnit",q+1,".txt",sep="")
    write.table(file=fnameA, PC1[,1:dims])
    write.table(file=fnameB, PCall[,1:dims])
    write.table(file=fnameC, PCu1[,1:dims])
    write.table(file=fnameD, PCuall[,1:dims])
    write.table(file=fnameE, PCuall.ctr[,1:dims])
    write.table(file=fnameF, PCuall.cor[,1:dims])
  }
  
}

out = list()
# PCA scores and loadings
out$scoreStart = PCu1[,c(1:dims)]
out$loadStart = PC1[,c(1:dims)]
out$rowpcoords =  PCuall[,c(1:dims)]  
out$colpcoords =  PCall[,c(1:dims)]

# PCA eigenvalues
sv = eg12$val/sqrt(nrows)
out$inertia_e=sv/(sum(sv))
out$sv = sv[c(1:dims)] 
out$levelnames = collabs
out$rownames = rowlabs
# Row contributions and correlations
out$rowctr=PCuall.ctr[,c(1:dims)]
#rownames(out$rowctr) = rowlabs
out$rowcor=PCuall.cor[,c(1:dims)]
#   rownames(out$rowcor) = rowlabs
out$nchunk = nchunk
out$disk = disk

if(disk==FALSE){
  out$allrowcoords=allCoordsU
  out$allcolcoords=allCoords
  out$allrowctr=allctrU
  out$allrowcor=allcorU
}
class(out)="i_pca"
return(out)
}
