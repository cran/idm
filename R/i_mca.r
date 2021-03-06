i_mca <- function(data1, data2=NULL,method=c("exact","live"), 
                  current_rank,nchunk = 2,ff = 0,disk=FALSE) {
  
  if(anyNA(data1)==TRUE | anyNA(data2)==TRUE ){stop("The data set should not contain missing data")}
  
  if(is.null(data2)==FALSE) {
    if(dim(data1)[2] != dim(data2)[2]){stop("The data sets must have the same number of columns/variables")}
  }
  
  if(is.null(data2)==TRUE){
    out = mjca(data1,lambda="indicator",ret=TRUE)
    outZ = transform_z(data1,is.weight=FALSE)
    out$m = nrow(data1)
    out$rowmass = outZ$r
    out$orgn = colMeans(outZ$SZ[1:nrow(data1),])
    
  }else{
    if(method=="exact"){
      out = h_exact_mca(data1 = data1, data2 = data2, current_rank=current_rank, nchunk=nchunk, disk=disk)  
    }
    
    if(method=="live"){
      out = r_live_mca(data1 = data1, data2 = data2, current_rank=current_rank, nchunk=nchunk, ff=ff,disk=disk)  
    }
  }
  class(out)="i_mca"
  return(out)
}
