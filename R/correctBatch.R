correctBatch <-
function(res11,designX,Ys,batchCL,ngrps,grps,Conf,nt,cGsub,comres,theme,type){
  XD=standardize(designX)
  Covars=cbind(rep(1,nrow(Ys)),XD)
  if(ngrps==2){resid=res11$u-tcrossprod(matrix(res11$coefficients[,2]),matrix(Covars[,2]))
  }else{resid=res11$u-tcrossprod(res11$coefficients[,2:ngrps],Covars[,2:ngrps]) }
  N=nrow(Ys); p=ncol(Ys); Nsamples=200
  mumut=matrix(colMeans(Ys),nrow=p,ncol=N,byrow=FALSE)
  psignal=t(res11$W%*%resid + mumut)
  set.seed(111)
  simD=lapply(1:Nsamples,function(x){psignal+matrix(rnorm(N*p,0,res11$Sig),N,p)})
  simData=apply(simplify2array(simD),c(1,2),median)
  ppccaDat=apply(simData,c(1,2),median)
  rownames(ppccaDat)=rownames(Ys)
  colnames(ppccaDat)=colnames(Ys)
  filname1=paste0(getwd(),"/correctBATCH/",Sys.Date(),"_ppccaCorrectedData.txt")
  write.table(ppccaDat,filname1,quote=TRUE,sep="\t",row.names=TRUE)
  ### pca on ppcca batch corrected data
  res4<-prcomp(ppccaDat, scale=TRUE )$x[,1:2]
  x=res4[,1]
  y=res4[,2]
  sdat5=data.frame(x=x,y=y,group=grps[batchCL])
  tdat5=data.frame(sdat5[,1:2],grps[batchCL],Conf)
  colnames(tdat5)<-c("PC-1","PC-2","Batch","Type")
  rownames(tdat5)<-rownames(Ys)
  filname4=paste0(getwd(),"/correctBATCH/",Sys.Date(),"_pca_After_PPCCA_correction.txt")
  write.table(tdat5,filname4,quote=TRUE,sep="\t",row.names=TRUE)
  if(nt==1){
    fig6<-ggplot(data= sdat5, aes(x, y,colour=cGsub[batchCL])) +
      geom_point(size=5.5,alpha=.7)+
      scale_shape_manual(guide=guide_legend(override.aes=aes(size=5)))+
      labs(title="",x="PC-1",y="PC-2") +theme_bw()+theme+
      scale_colour_manual("Batch",labels=grps,values=cGsub)+
      guides(colour = guide_legend(override.aes = list(linetype=0)))
  }else{
    if(nt>8){
      z=(Conf-min(Conf))/(max(Conf)-min(Conf))
      z=((z+.15)*10)
      fig6<-ggplot(data= sdat5, aes(x, y,colour=cGsub[batchCL])) +
        geom_point(size=z,alpha=.7)+
        scale_shape_manual(values=fch,guide=guide_legend(override.aes=aes(size=z)))+
        labs(title="",x="PC-1",y="PC-2") +theme_bw()+theme+
        scale_colour_manual("Batch",labels=grps,values=cGsub)+
        guides(colour = guide_legend(override.aes = list(linetype=0)))
    }else{
      pch=c(15:22)
      fch=pch[1:nt]
      lb=paste0('bio-',type)
      fig6<-ggplot(data= sdat5, aes(x, y,colour=cGsub[batchCL])) +
        geom_point(aes(shape=factor(fch[Conf],labels=lb)),size=5.5,alpha=.7)+
        scale_shape_manual(values=fch,guide=guide_legend(override.aes=aes(size=5)))+
        labs(title="",x="PC-1",y="PC-2") +theme_bw()+theme+
        scale_colour_manual("Batch",labels=grps,values=cGsub)+
        guides(colour = guide_legend(override.aes = list(linetype=0)))
    }
  }
  filname5=paste0(getwd(),"/correctBATCH/",Sys.Date(),"_pca_After_PPCCA_correction.png") 
  suppressMessages(ggsave(filename=filname5,plot = fig6))
  ############## ppcca vs combat correction
  x=c(t(comres))  # check if you need scaled data instead? No, GaussianSpray allows genes to have var!=1
  y=c(ppccaDat)  
  sdat7=tdat7<-data.frame(x=x,y=y)
  colnames(tdat7)<-c("CombatPredictedData","PPCCApredictedData")
  filname6=paste0(getwd(),"/correctBATCH/",Sys.Date(),"_ComBat_vs_PPCCA_predictedData.txt")
  write.table(tdat7,filname6,quote=TRUE,sep="\t",row.names=TRUE)
  fig8<-ggplot(data=sdat7, aes(x, y,colour=c("orange")[rep(1,length(y))])) + 
    geom_point(aes(shape=factor(16)),size=3,alpha=.2)+scale_colour_manual("",values="orange")+theme(legend.key = element_blank())+
    labs(title=paste0("Correlation: ",round(cor(x,y),2)),x="Combat-predicted data",y="PPCCA-predicted data") +theme_bw()+theme+
    theme(legend.position = c(10000, 10000))+theme(legend.key = element_blank())+theme(legend.title=element_blank())+
    guides(colour = guide_legend(override.aes = list(linetype=0)))+
    coord_fixed(ratio = 1.15)
  filname7=paste0(getwd(),"/correctBATCH/",Sys.Date(),"_ComBat_vs_PPCCA_predictedData.png")
  suppressMessages(ggsave(filename=filname7,plot = fig8))
  rm(x,y,psignal,simD,res4,sdat7,tdat7,fig8,comres)
  gc()
  ppccaDat}
