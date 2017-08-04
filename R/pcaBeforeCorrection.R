pcaBeforeCorrection <- function(rs1,grps,cGsub,batchCL,Conf,type,Ys,theme){
  #### proportion of variation
  variation <- c(rs1$sdev^2/sum(rs1$sdev^2))
  x <- 1:length(variation)
  y <- variation
  variat <- data.frame(x=x, y=y)
  colnames(variat)<-c("numberPCs","PropVariation")
  rownames(variat)<-paste0("PC-",1:length(variation))
  write.table(variat,paste0("pcaBeforeCorrection/",Sys.Date(),"_ProportionVariationPCAbeforeCorrection.txt"),quote=TRUE,sep="\t",row.names=TRUE)
  prop <- ggplot(data=variat, aes(x, y)) + geom_point(size=3)+geom_line()+
    labs(title=" ",x="Number of PCs",y="Proportion of variation")+theme_bw()+theme
  suppressMessages(ggsave(filename ="pcaBeforeCorrection/ProportionVariationPCAbeforeCorrection.pdf", device=cairo_pdf,plot = prop))
  #### pca scores
  x <- rs1$x[,1]
  y <- rs1$x[,2]
  sdat1 <- data.frame(x=x,y=y,group=grps[batchCL])
  tdat1<-data.frame(sdat1[,1:2],grps[batchCL],Conf)
  colnames(tdat1)<-c("PC-1","PC-2","Batch","Type")
  rownames(tdat1)<-rownames(Ys)
  write.table(tdat1,paste0("pcaBeforeCorrection/",Sys.Date(),"_pcaBeforeCorrection.txt"),quote=TRUE,sep="\t",row.names=TRUE)
  
  if(length(type)==1){
    fig1<-ggplot(data= sdat1, aes(x, y,colour=cGsub[batchCL])) +
      geom_point(size=5.5,alpha=.85)+
      scale_shape_manual(guide=guide_legend(override.aes=aes(size=5)))+
      labs(title="",x="PC-1",y="PC-2") +theme_bw()+theme+
      scale_colour_manual("Batch",labels=grps,values=cGsub)+
      guides(colour = guide_legend(override.aes = list(linetype=0)))
  }else{
    if(length(type)>8){
      z=(Conf-min(Conf))/(max(Conf)-min(Conf))
      z=((z+.15)*10)
      fig1<-ggplot(data= sdat1, aes(x, y,colour=cGsub[batchCL])) +
        geom_point(size=z,alpha=.5)+
        scale_shape_manual(values=fch,guide=guide_legend(override.aes=aes(size=z)))+
        labs(title="",x="PC-1",y="PC-2") +theme_bw()+theme+
        scale_colour_manual("Batch",labels=grps,values=cGsub)+
        guides(colour = guide_legend(override.aes = list(linetype=0)))
    }else{
      pch=c(15:22)
      fch=pch[1:length(type)]
      lb=paste0('bio-',type)
      fig1<-ggplot(data= sdat1, aes(x, y,colour=cGsub[batchCL])) +
        geom_point(aes(shape=factor(fch[Conf],labels=lb)),size=5.5,alpha=.7)+
        scale_shape_manual(values=fch,guide=guide_legend(override.aes=aes(size=5)))+
        labs(title="",x="PC-1",y="PC-2") +theme_bw()+theme+
        scale_colour_manual("Batch",labels=grps,values=cGsub)+
        guides(colour = guide_legend(override.aes = list(linetype=0)))
    }
  }
  suppressMessages(ggsave(filename ="pcaBeforeCorrection/pcaBeforeCorrection.png",plot = fig1))
  #### pairs plot for pca
  pd<-rs1$x
  colnames(pd)<-paste0("PC-",1:ncol(rs1$x))
  write.table(pd,paste0("pcaBeforeCorrection/",Sys.Date(),"pairspcaBeforeCorrection.txt"),quote=TRUE,sep="\t",row.names=TRUE)
  pdf(file="pcaBeforeCorrection/pairsplotpcaBeforeCorrection.pdf",onefile=TRUE)
  pairs(pd[,1:10],col=cGsub[batchCL],pch=15)
  dev.off()
  rm(pd,fig1,tdat1,sdat1,prop,variat,variation)
  gc()}
