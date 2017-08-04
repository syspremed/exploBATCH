correctComBat <-
function(Ys,batchCL,grps,cGsub,Conf,theme,type){
  comdat<-t(Ys); 
  comres<-ComBat(dat=comdat,batch=batchCL,mod=NULL, par.prior=TRUE, prior.plots=FALSE) # run combat
  write.table(t(comres),paste0("correctComBat/",Sys.Date(),"combatCorrectedData.txt"),quote=TRUE,sep="\t",row.names=TRUE)
  ### pca on combat batch corrected
  res5<-prcomp(t(comres),scale=TRUE)
  x <- res5$x[,1]
  y <- res5$x[,2]
  sdat6 <- data.frame(x=x,y=y,group=grps[batchCL])
  tdat6<-data.frame(sdat6[,1:2],grps[batchCL],Conf)
  colnames(tdat6)<-c("PC-1","PC-2","Batch","Type")
  rownames(tdat6)<-rownames(Ys)
  write.table(tdat6,paste0("correctComBat/",Sys.Date(),"pca_After_Combat_Correction.txt"),quote=TRUE,sep="\t",row.names=TRUE)
  if(length(type)==1){
    fig7<-ggplot(data= sdat6, aes(x, y,colour=cGsub[batchCL])) +
      geom_point(size=5.5,alpha=.7)+
      scale_shape_manual(guide=guide_legend(override.aes=aes(size=5)))+
      labs(title="",x="PC-1",y="PC-2") +theme_bw()+theme+
      scale_colour_manual("Batch",labels=grps,values=cGsub)+
      guides(colour = guide_legend(override.aes = list(linetype=0)))
  }else{
    if(length(type)>8){
      z=(Conf-min(Conf))/(max(Conf)-min(Conf))
      z=((z+.15)*10)
      fig7<-ggplot(data= sdat6, aes(x, y,colour=cGsub[batchCL])) +
        geom_point(size=z,alpha=.7)+
        scale_shape_manual(values=fch,guide=guide_legend(override.aes=aes(size=z)))+
        labs(title="",x="PC-1",y="PC-2") +theme_bw()+theme+
        scale_colour_manual("Batch",labels=grps,values=cGsub)+
        guides(colour = guide_legend(override.aes = list(linetype=0)))
    }else{
      pch=c(15:22)
      fch=pch[1:length(type)]
      lb=paste0('bio-',type)
      fig7<-ggplot(data= sdat6, aes(x, y,colour=cGsub[batchCL])) +
        geom_point(aes(shape=factor(fch[Conf],labels=lb)),size=5.5,alpha=.7)+
        scale_shape_manual(values=fch,guide=guide_legend(override.aes=aes(size=5)))+
        labs(title="",x="PC-1",y="PC-2") +theme_bw()+theme+
        scale_colour_manual("Batch",labels=grps,values=cGsub)+
        guides(colour = guide_legend(override.aes = list(linetype=0)))
    }
  }
  suppressMessages(ggsave(filename ="correctComBat/pca_After_Combat_Correction.png",plot = fig7))
  comres}
