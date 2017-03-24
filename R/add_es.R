add_es <- function(eg,eg2,current_rank,ff=0,method=c("esm","isvd")){
  if(method=="esm"){
    if (missing("current_rank")) {
      #current_rank = full rank
      current_rank =  length(eg$d)
    }
 
    out = add_eig(eg, eg2, current_rank)}
  else{
    if (is.null(eg$m)) {
      m = dim(eg$u)[1]
    } 
    else
    {
      m = eg$m
    }
    B = eg2
    if (missing("current_rank")) {
      #current_rank = full rank
      current_rank =  dim(B)[2]
    }
    out = add_svd(eg,B,m,current_rank,ff = ff)
  }
  return(out)
}