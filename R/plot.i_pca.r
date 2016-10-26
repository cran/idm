plot.i_pca<-function(x,dims=c(1,2),what=c(TRUE,TRUE),dataname=NULL,
                     labels=NULL, animation=TRUE,frames=10,zoom=TRUE,movie_format="gif",...){
  
  obj <- x
  
  nch =length(obj$allcolcoord)
  
  if(nch==0){
    obj$allcolcoord=list()
    obj$allcolcoord[[1]]=obj$colpcoord  
    nch=1
  }
  
  for(j in 1:nch ){
    dist2var=apply(obj$allcolcoord[[j]],1,function(x)(sqrt(sum(x^2))))
    obj$allcolcoord[[j]]=obj$allcolcoord[[j]]/dist2var
  }
  
  #default dataname
  if (is.null(dataname)){
    dataname = deparse(substitute(x))
  }
  
  #default labels
  if (is.null(labels)){
    labels = obj$levelnames
  }
  
  if(animation==TRUE){
    outmo=all_frame_make(obj=obj,dims=dims,nfrs=frames,is.PCA=TRUE)
    
    #attributes
    if (what[2] == TRUE) {
      
      movieNameA=paste("iPCA_",dataname,"atts_","movie",sep="")#,movie_format
      ani_plot(moname=movieNameA,outmo,nfrs=frames,labs=labels,att=TRUE,pca=TRUE,zoom=zoom,
               movie_format=movie_format)
    }
    
    #objects
    if (what[1] == TRUE) {
      movieNameO=paste("iPCA_",dataname,"obs_","movie",sep="")#movie_format,
      ani_plot(moname=movieNameO,outmo,nfrs=frames,labs=labels,att=FALSE,pca=TRUE,zoom=zoom,
               movie_format=movie_format)
    }
  }else{
    stpl=static_plot(obj,dims=dims,what=what,labs=labels,pca=TRUE)
    return(stpl)
  }
  
  
}
