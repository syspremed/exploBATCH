assessComBat <-
function(rerun1,qopt,nt,rerun12,theme){
  if(nt!=1){
    #### Detect biological effect
    pnames=colnames(rerun12$coefficients)
    for(i in 2:ncol(rerun12$coefficients)){ 
      blci<-data.frame(matrix(paste0("pPC-",1:qopt)),matrix(rerun12$coefficients[,i]),rerun12$coeffCI[,(2*i-1):(2*i)])
      colnames(blci)<-c("x","y","ylo","yhi")
      bldat<-blci[,-1]
      rownames(bldat)<-blci[,1]
      colnames(bldat)<-c("Effect","LowerCI","UpperCI")
      filnamec1=paste0(getwd(),"/assessComBat/",Sys.Date(),"_BioEFFECT.txt")
      write.table(bldat,filnamec1,quote=TRUE,sep="\t",row.names=TRUE)
      bpl <- ggplot(blci, aes(x=x, y=y, ymin=ylo, ymax=yhi))+geom_pointrange()+labs(title="",x="",y="Biological effect")+
        theme_bw()+theme+coord_flip()+geom_hline(aes(yintercept=0),lty=5,col=2,size=1.5)
      suppressMessages(ggsave(filename=paste0(getwd(),"/assessComBat/BioEFFECT-",pnames[i],".png"),plot = bpl))}}
  #### Detect batch effect
  lci<-data.frame(matrix(paste0("pPC-",1:qopt)),matrix(rerun1$coefficients[,2]),rerun1$coeffCI[,3:4])
  colnames(lci)<-c("x","y","ylo","yhi")
  ldat<-lci[,-1]
  rownames(ldat)<-lci[,1]
  colnames(ldat)<-c("Effect","LowerCI","UpperCI")
  write.table(ldat,paste0("assessComBat/",Sys.Date(),"batchEffect.txt"),quote=TRUE,sep="\t",row.names=TRUE)
  pl<-ggplot(lci, aes(x=x, y=y, ymin=ylo, ymax=yhi))+geom_pointrange()+labs(title="",x="",y="Batch effect")+
    theme_bw()+theme+coord_flip()+geom_hline(aes(yintercept=0),lty=5,col=2,size=1.5)
  suppressMessages(ggsave(filename="assessComBat/BatchEffect.png",plot = pl))
}
