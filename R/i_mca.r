i_mca <- function(data1, data2,method=c("exact","live"), nchunk = 2,f = 0,disk=FALSE) {
  
  #  source("h_exact_mca.r")
  #  source("r_live_mca.r")
  
  if(method=="exact"){
    out=h_exact_mca(data1=data1, data2=data2,nchunk=nchunk,disk=disk)  
  }
  
  if(method=="live"){
    out=r_live_mca(data1=data1, data2=data2,nchunk=nchunk,f=f,disk=disk)  
  }
  class(out)="i_mca"
  return(out)
}