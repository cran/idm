i_pca <- function(data1, data2 = NULL, current_rank, nchunk = 2, disk=FALSE) {
  ## data1 the starting matrix; data2 the upcoming data
  ## nchunk: number of chunks to split data2
  ## disk: store output to use for graphics
  if(anyNA(data1)==TRUE | anyNA(data2)==TRUE ){stop("The data set should not contain missing data")}
  ## This is equivalent to PCA on the covariance matrix 
  ## (i.e. on mean centered data)
  
  if(is.null(data2)==FALSE) {
    if(dim(data1)[2] != dim(data2)[2]){stop("The data sets must have the same number of columns/variables")}
  }
  
  if (missing("current_rank")) {
    #current_rank = full rank
    current_rank = dim(data1)[2]
  }
  
  tdata2 = FALSE
  if(is.null(data2)==TRUE){
    tdata2 = TRUE #keep temp state of data2
    n=nrow(data1)
    n1s=1:ceiling(n/2)
    
    data2=data1[-n1s,]
    data1=data1[n1s,]
  }
  if(disk==TRUE){
    suppressWarnings(dir.create("./PCstory"))
  }else{
    allCoords=list()
    allCoordsU=list()
    allctr=list()
    allcor=list()
    allctrU=list()
    allcorU=list()
  }
  
  #rowlabs = c(rownames(data1),rownames(data2))
  collabs = colnames(data1)
  nrows = nrow(rbind(data1,data2))
  ncols = ncol(data1)
  nrows1 = nrow(data1)
  
  dims = current_rank #=ncols
  
  if ((length(nchunk) > 1 ) & (sum(nchunk) != nrow(data2))) {
    stop("\nchunk blocks do not match the size of 'data2'")
  }
  n1=nrow(data1)
  data1=scale(data1,center=F,scale=F)#*sqrt(n1/(n1-1))
  #print(head(data1[1:6,1:5]))
  ########################################################################
  ########################################################################
  #Try that everything works without scaling
  ########################################################################
  
  eg1 = do_es(data1)
  
  PC1 = (1/sqrt(n1))*eg1$v[,1:current_rank]%*%diag(eg1$d[1:current_rank])
  PCu1 = eg1$u[,1:current_rank]%*%diag(eg1$d[1:current_rank]) 

  ## insert ctr comps
  
  eig=(eg1$d[1:current_rank] * (1/sqrt(n1)))^2
  #print(eig)
  
  #print(scale(data1,center=T,scale=F)[1:6,1:3])
  #dist2.ind <- as.vector(tcrossprod(as.matrix(data1^2 ), t(rep(1, ncol(data1)))))
  cen_data1=scale(data1,center=T,scale=F)
  dist2.ind <- apply(cen_data1^2,1,sum)
  sum_sq_d1 = apply(data1^2,2,sum)
  m1=apply(data1,2,mean)
  dist2.var = (1/n1)*(sum_sq_d1-nrow(data1)* m1^2)
  #dist2.var <- (1/n1)*apply(cen_data1^2,2,sum)
  #print(mydist2.var)
  #print(dist2.var)
  # mydist2.var=as.vector(crossprod(rep(1, nrow(data1)), as.matrix(cen_data1^2 * (1/n1))))
  
  PC1.ctr = (PC1^2)%*% pseudoinverse(diag(1/eig))
  
  PC1.cor= PC1/sqrt(dist2.var)
  #print(PC1.cor[1:6,1:2])
  
  PCu1.cor <- PCu1^2/dist2.ind
  #print(PCu1.cor[1:6,1:2])
  PCu1.ctr <- t(t(PCu1^2 * (1/n1))/eig)
  # print(PCu1.ctr[1:6,1:2])
  ########################################################################
  ########################################################################
  ########################################################################
  
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
    allctr[[1]] = PC1.ctr[,c(1:dims)]
    allcor[[1]] = PC1.cor[,c(1:dims)]
  }
  
  if (length(nchunk) == 1)
  {
    if(floor(nrow(data2)/nchunk) < 2){stop("There are blocks with zero rows. Please reduce the number of chunks.")}
  } else {
    if(any(nchunk < 1) == TRUE){stop("Block size should be greater than 1.")}
  }
  out.split = mat_split(data2, (nchunk)) 
  mat.story = out.split$splitMat
  
  #if block sizes are given, switch back to number
  if (length(nchunk) > 1) {
    nchunk = length(nchunk)
  } 
  
  for (q in 1:length(mat.story)) {
    
    mat.chu = data.matrix(mat.story[[q]])
    nchu=nrow(mat.chu)
    #mat.chu=scale(mat.chu,center=F,scale=T)*sqrt(nchu/(nchu-1))
    ### coordinate computation
    ## column computation (modalities)
    #     print(q)
    #     if (q > 1) {
    #       #######################
    #       PCu1 = PCuall
    #       PC1 = PCall
    #       #######################
    #     }
    
    ########################################################################
    ########################################################################
    ########################################################################
    #   if (dim(mat.chu)[1] < dim(mat.chu)[2]) {
    #      mat.chu = t(mat.chu)
    #    }
    nchu = nrow(mat.chu)
    
    mat.chu = scale(mat.chu,center=F,scale=F)#*sqrt(nchu/(nchu-1))
    sum_sq_dchu = apply(mat.chu^2,2,sum)
    eg2 = do_es(mat.chu)
    
    eg12 = add_es(eg1, eg2, method="esm", current_rank)
    n12=eg12$m
    m12=eg12$orgn
    
    n2=(nrow(mat.chu))
    
    PCall = (1/sqrt(n12))* eg12$v%*%diag(eg12$d) 
    PCuall = eg12$u%*%diag(eg12$d) 
    # print(nrow(PCuall))
    
    eig=(eg12$d * (1/sqrt(n12)))^2
    # print(eig)
    dist2_12.ind = apply(PCuall^2,1,sum)
    
    #dist2chu.var <- (1/n2)*apply(cen_mat.chu^2,2,sum)
    sum_sq_d12=sum_sq_d1+sum_sq_dchu
    dist2_12.var = as.vector((1/n12)*((sum_sq_d12)-n12* m12^2))
    PCall.ctr = t(t(PCall^2)/eig)*100
    PCall.cor= (PCall / sqrt(dist2_12.var))^2
    
    PCuall.cor <- (PCuall^2)/dist2_12.ind
    PCuall.ctr <- t(t(PCuall^2 * (1/n12))/eig)*100
    
    # PCuall.ctr <- t(t(coord.ind^2 * row.w/sum(row.w))%*% diag(eig))
    PCall = sign_match(PC1, PCall)
    PCuall = sign_match(PCu1, PCuall)
    ########################################################################
    ########################################################################
    ########################################################################
    if(disk==FALSE){      
      allCoords[[q+1]]=PCall[,c(1:dims)]
      allCoordsU[[q+1]]=PCuall[,c(1:dims)]
      allctrU[[q+1]] = PCuall.ctr[,c(1:dims)]
      allcorU[[q+1]] = PCuall.cor[,c(1:dims)]
      allctr[[q+1]] = PCall.ctr[,c(1:dims)]
      allcor[[q+1]] = PCall.cor[,c(1:dims)]
    } 
    
    eg1 = eg12
    dist2.ind=dist2_12.ind
    #dist2.var=dist2_12.var
    sum_sq_d1=sum_sq_d12
    n1=n12
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
  out$rowpcoord =  PCuall[,c(1:dims)]  
  out$colpcoord =  PCall[,c(1:dims)]
  out$eg=eg12
  # PCA eigenvalues
  sv = eg12$d/sqrt(nrows)
  out$sv = sv[c(1:dims)]
  if (current_rank == ncols) {
    out$inertia.e= eg12$d^2/(sum(eg12$d^2))
  } else {
    out$inertia.e= sv^2/ncols 
  }
  out$levelnames = collabs
  # out$rownames = rowlabs
  # Row contributions and correlations
  out$rowctr=PCuall.ctr[,c(1:dims)]
  out$colctr=PCall.ctr[,c(1:dims)]
  #rownames(out$rowctr) = rowlabs
  out$rowcor=PCuall.cor[,c(1:dims)]
  out$colcor=PCall.cor[,c(1:dims)]
  out$nchunk = nchunk
  out$disk = disk
  
  if((disk==FALSE) & (tdata2==FALSE)) {
    out$allrowcoord=allCoordsU
    out$allcolcoord=allCoords
    out$allrowctr=allctrU
    out$allcolctr=allctr
    out$allrowcor=allcorU
    out$allcolcor=allcor
  }
  class(out)="i_pca"
  return(out)
}
