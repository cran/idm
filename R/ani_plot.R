ani_plot <- function(outmo,nfrs=25,moname="mymovie",labs,att=TRUE,pca=TRUE,contrib,zoom,movie_format="gif",binary){
  
  oopt <- animation::ani.options(interval = 0.1)
  if(att==TRUE){
    if(zoom==FALSE){
      xr=range(unlist(lapply(1:outmo$nchunk,function(jj) {outmo[[jj]]$xr })))
      yr=range(unlist(lapply(1:outmo$nchunk,function(jj) {outmo[[jj]]$yr })))
    }else{xyr=active_zoom(outmo,nframes=nfrs,att=T)}
  }else{
    if(zoom==FALSE){
      xr=range(unlist(lapply(1:outmo$nchunk,function(jj) {outmo[[jj]]$uxr })))
      yr=range(unlist(lapply(1:outmo$nchunk,function(jj) {outmo[[jj]]$uyr })))  
    }else{xyr=active_zoom(outmo,nframes=nfrs,att=F)}
  }
  
  FUN2 <- function(nfrs,nchunk,xr,yr,contrib) {
    for (chu in 1:nchunk){
      if(zoom==FALSE){
        lapply(seq(1,nfrs, by = 1), function(i) {
          plot_fun(outmo,chu,i,xr,yr,lab=labs,att=att,pca=pca,contrib=contrib,binary)
          animation::ani.pause()})
      }else{
        lapply(seq(1,nfrs, by = 1), function(i) {
          plot_fun(outmo,chu,i,xyr[[chu]]$xr[i,],xyr[[chu]]$yr[i,],lab=labs,att=att,pca=pca,contrib=contrib,binary)
          animation::ani.pause()})
        
      }
    }
  } 
  
  if (file.exists(moname)) file.remove(moname)
  if(movie_format=="gif"){
    saveGIF(FUN2(nfrs,outmo$nchunk,xr,yr,contrib), interval = 0.1,movie.name=paste(moname,".",movie_format,sep=""))
  }else{
    if(att==TRUE){
      frame_name="att_frame"
    }else{
      frame_name="obs_frame"  
    }
    nmax = nfrs*outmo$nchunk
    saveLatex(FUN2(nfrs,outmo$nchunk,xr,yr,contrib), img.name = frame_name, ani.opts = "controls,width=0.95\\textwidth",
              latex.filename = ifelse(interactive(), paste(moname,".tex",sep=""), ""),
              interval = 0.1, nmax = nmax, ani.dev = "pdf", ani.type = "pdf", ani.width = 7,
              ani.height = 7,documentclass = paste("\\documentclass{article}",
                                                    "\\usepackage[papersize={7in,7in},margin=0.3in]{geometry}",
                                                    sep = "\n"))
  }
}


