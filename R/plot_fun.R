plot_fun<-function(outmo,chu,i,xr,yr,lab,att=TRUE,pca=TRUE,contrib="none",binary){
  #  require("ggplot2")
  #set to NULL to avoid compiler NOTE 
  x <- NULL
  y <- NULL
  ctr <- NULL
  p=ncol(outmo[[1]]$dat$Xs)
  # suppressWarnings(if(lab=="foo"){lab=paste("v_",1:p,sep="")})
  if(att==TRUE){
    if(pca==TRUE){
      cdf=circle_fun()
      df=data.frame(x=outmo[[chu]]$dat$Xs[i,],y=outmo[[chu]]$dat$Ys[i,],
                    lab=lab)
      #,ctr=outmo[[chu]]$bfcolctr[i,],
      #cor=outmo[[chu]]$bfcolcor[i,])
      a=ggplot(data=df,aes(x=x,y=y))#+xlim(xr)+ylim(yr)
      a=a+geom_text(data=df,aes(label=lab))+ xlab("")+ylab("") #,size=ctr
      a=a+geom_point(data=cdf,aes(x=x,y=y),size=.05)
      a=a+geom_segment(data=df,aes(x=0,xend=x,y=0,yend=y),size=.5,alpha=.5)
      a=a+geom_vline(xintercept=0)+geom_hline(yintercept=0)+theme_bw()
    }else{
      #when matrix is presence/absence
      if (binary == TRUE) {
        df=data.frame(x=outmo[[chu]]$dat$Xs[i,seq(2,p,2)],y=outmo[[chu]]$dat$Ys[i,seq(2,p,2)],ctr=outmo[[chu]]$bfcolctr[i,seq(2,p,2)],cor=outmo[[chu]]$bfcolcor[i,seq(2,p,2)],lab=lab)
      } else {
      df=data.frame(x=outmo[[chu]]$dat$Xs[i,],y=outmo[[chu]]$dat$Ys[i,],ctr=outmo[[chu]]$bfcolctr[i,],
                    cor=outmo[[chu]]$bfcolcor[i,],lab=lab)
      }
      a=ggplot(data=df,aes(x=x,y=y))+xlim(xr)+ylim(yr)
      a=a+geom_vline(xintercept=0)+geom_hline(yintercept=0)+theme_bw()
      a=a+theme(legend.position="none")
      if (contrib == "none") {
        a=a+geom_text(data=df,aes(label=lab))+ xlab("")+ylab("")
      }
      if (contrib == "cor")
      {
        a=a+geom_text(data=df,aes(label=lab,size=cor))+ xlab("")+ylab("")
      }
      
      if (contrib == "ctr") {
        a=a+geom_text(data=df,aes(label=lab,size=ctr))+ xlab("")+ylab("")
      }
    }
  }else{
    df=data.frame(x=outmo[[chu]]$udat$Xs[i,],y=outmo[[chu]]$udat$Ys[i,])
    a=ggplot(data=df,aes(x=x,y=y))+xlim(xr)+ylim(yr)
    a=a+geom_point(size=.5)+ xlab("")+ylab("")
    a=a+geom_vline(xintercept=0)+geom_hline(yintercept=0)+theme_bw()
  }
  
  #  return(a)
  print(a)
  
}
