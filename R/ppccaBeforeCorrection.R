ppccaBeforeCorrection <-
function(res1,grps,cGsub,batchCL,Conf,type,Ys,theme){
  #### number of ppcs
  x=c(2:(length(res1$BIC)+1))
  y=res1$BIC
  bic <- bicp <- data.frame(x=x,y=y)
  colnames(bicp)<-c("Number-pPCs","BIC")
  rownames(bicp)<-paste0("pPC-",2:(length(res1$BIC)+1))
  write.table(bicp,paste0("ppccaBeforeCorrection/",Sys.Date(),"_NumberOfpPCs.txt"),quote=TRUE,sep="\t",row.names=TRUE)
  fig2<-ggplot(data=bic, aes(x, y)) + geom_point(size=3)+geom_line()+labs(title=" ",x="Number of pPCs", 
        y="BIC")+theme_bw()+theme+geom_vline(xintercept=res1$q,lty=5,col=2,size=1.5)
  suppressMessages(ggsave(filename ="ppccaBeforeCorrection/NumberOfpPCs.png",,plot=fig2))
  #### ppcca scoresplot before batch correction 
  x=res1$u[1,]
  y=res1$u[2,]
  sdat2=data.frame(x=x,y=y,group=grps[batchCL])
  tdat2=data.frame(sdat2[,1:2],grps[batchCL],Conf)
  colnames(tdat2)<-c("pPC-1","pPC-2","Batch","Type")
  rownames(tdat2)<-rownames(Ys) 
write.table(tdat2,paste0("ppccaBeforeCorrection/",Sys.Date(),"ppccaBeforeCorrection.txt"),quote=TRUE,sep="\t",row.names=TRUE)
  if(length(type)==1){
    fig3<-ggplot(data= sdat2, aes(x, y,colour=cGsub[batchCL])) +
      geom_point(size=5.5,alpha=.7)+
      scale_shape_manual(guide=guide_legend(override.aes=aes(size=5)))+
      labs(title="",x="pPC-1",y="pPC-2") +theme_bw()+theme+
      scale_colour_manual("Batch",labels=grps,values=cGsub)+
      guides(colour = guide_legend(override.aes = list(linetype=0)))
  }else{
    if(length(type)>8){
      z=(Conf-min(Conf))/(max(Conf)-min(Conf))
      z=((z+.15)*10)
      fig3<-ggplot(data= sdat2, aes(x, y,colour=cGsub[batchCL])) +
        geom_point(size=z,alpha=.7)+
        scale_shape_manual(values=fch,guide=guide_legend(override.aes=aes(size=z)))+
        labs(title="",x="pPC-1",y="pPC-2") +theme_bw()+theme+
        scale_colour_manual("Batch",labels=grps,values=cGsub)+
        guides(colour = guide_legend(override.aes = list(linetype=0)))
    }else{
      pch=c(15:22)
      fch=pch[1:length(type)]
      lb=paste0('bio-',type)
      fig3<-ggplot(data= sdat2, aes(x, y,colour=cGsub[batchCL])) +
        geom_point(aes(shape=factor(fch[Conf],labels=lb)),size=5.5,alpha=.7)+
        scale_shape_manual(values=fch,guide=guide_legend(override.aes=aes(size=5)))+
        labs(title="",x="pPC-1",y="pPC-2") +theme_bw()+theme+
        scale_colour_manual("Batch",labels=grps,values=cGsub)+
        guides(colour = guide_legend(override.aes = list(linetype=0)))
    }
  }
  suppressMessages(ggsave(filename ="ppccaBeforeCorrection/ppccaBeforeCorrection.png",plot = fig3))
  #### pairs plot for ppcca
  pdat<-t(res1$u)
  colnames(pdat)<-paste0("pPC-",1:res1$q)
  write.table(pdat,paste0("ppccaBeforeCorrection/",Sys.Date(),"_pairsppccaBeforeCorrection.txt"),quote=TRUE,sep="\t",row.names=TRUE)
  pdf(file="ppccaBeforeCorrection/pairsplotppccaBeforeCorrection.pdf",onefile=TRUE)
  pairs(pdat[,1:res1$q],col= cGsub[batchCL],pch=15)
  dev.off()
  rm(pdat,fig3,sdat2,tdat2,fig2,bicp,bic)
  gc()
  }
