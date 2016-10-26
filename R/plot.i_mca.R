plot.i_mca<-function(x,dims=c(1,2),what=c(TRUE,TRUE),contrib="none",dataname=NULL,
                     labels=NULL,animation=TRUE,frames=10,zoom=TRUE,movie_format="gif",binary=FALSE, ...){
  
  # require(animation)
  #  source("static_plot.R")
  #  source("all_frame_make.R")
  #  source("ani_plot.R")
  obj <- x
  #default dataname
  if (is.null(dataname)){
    dataname = deparse(substitute(x))
  }
  
  #default labels
  if (is.null(labels)){
    labels =  obj$levelnames
    #case of presence/absence dataset
    if (binary == TRUE) {
      labels = obj$levelnames[seq(2,length(obj$levelnames),2)]
      #remove last two characters
      labels = substr(labels, 1, nchar(labels)-2)
    }
  }
  if(animation==TRUE){
    outmo=all_frame_make(obj=obj,dims=dims,nfrs=frames,is.PCA=FALSE)  
    
    #attributes
    if (what[2] == TRUE) {
      movieNameA=paste("iMCA_",dataname,"atts_","movie",sep="")
      ani_plot(moname=movieNameA,outmo,nfrs=frames,labs=labels,att=TRUE,pca=FALSE,contrib=contrib,zoom=zoom,
               movie_format=movie_format,binary)
    }
    
    #objects
    if (what[1] == TRUE) {
      movieNameO=paste("iMCA_",dataname,"obs_","movie",sep="")
      ani_plot(moname=movieNameO,outmo,nfrs=frames,labs=labels,att=FALSE,pca=FALSE,contrib=contrib,zoom=zoom,
               movie_format=movie_format,binary)
    }
  }else{
    stpl=static_plot(obj, dims=dims, what=what,labs=labels,pca=FALSE,contrib=contrib,binary=binary)
    return(stpl)
  }
}

