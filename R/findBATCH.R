findBATCH <-
function(res2,qopt,nt,rerun,theme){
  #### Detect biological effect
  if(nt!=1){
    pnames=colnames(rerun$coefficients)
    for(i in 2:ncol(rerun$coefficients)){ 
      blci<-data.frame(matrix(paste0("pPC-",1:qopt)),matrix(rerun$coefficients[,i]),rerun$coeffCI[,(2*i-1):(2*i)])
      colnames(blci)<-c("x","y","ylo","yhi")
      bldat<-blci[,-1]
      rownames(bldat)<-blci[,1]
      colnames(bldat)<-c("Effect","LowerCI","UpperCI")
      filnamef1=paste0(getwd(),"/findBATCH/",Sys.Date(),"_BioEFFECT.txt")
      write.table(bldat,filnamef1,quote=TRUE,sep="\t",row.names=TRUE)
      bpl <- ggplot(blci, aes(x=x, y=y, ymin=ylo, ymax=yhi))+geom_pointrange()+labs(title="",x="",y="Biological effect")+
        theme_bw()+theme+coord_flip()+geom_hline(aes(yintercept=0),lty=5,col=2,size=1.5)
      suppressMessages(ggsave(filename=paste0(getwd(),"/findBATCH/BioEFFECT-",pnames[i],".png"),plot = bpl))}}
  #### Detect batch effect
  lci<-data.frame(matrix(paste0("pPC-",1:qopt)),matrix(res2$coefficients[,2]),res2$coeffCI[,3:4])
  LCI<-matrix((lci[,3]*lci[,4])>0)
  colnames(lci)<-c("x","y","ylo","yhi")
  ldat<-lci[,-1]
  rownames(ldat)<-lci[,1]
  colnames(ldat)<-c("Effect","LowerCI","UpperCI")
  write.table(ldat,paste0("findBATCH/",Sys.Date(),"batchEffect.txt"),quote=TRUE,sep="\t",row.names=TRUE)
  pl <- ggplot(lci, aes(x=x, y=y, ymin=ylo, ymax=yhi))+geom_pointrange()+labs(title="",x="",y="Batch effect")+
    theme_bw()+theme+coord_flip()+geom_hline(aes(yintercept=0),lty=5,col=2,size=1.5)
  suppressMessages(ggsave(filename="findBATCH/BatchEffect.png",plot = pl))
  LCI}
