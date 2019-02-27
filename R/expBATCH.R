expBATCH <-function(D,batchCL,Conf=NA,mindim=2,maxdim=3,method="ppcca",scale="unit",SDselect=0){
  if(missing(D)){stop("Expression data are required.\n")}
  if(missing(batchCL)){stop("At least a batch variable is required.\n ")}
  if(nrow(D)!= length(batchCL)){stop("Expression data and batch variable should be matching in samples.\n")}
  if((mindim<2)|(maxdim<2)){stop("Minimum number of PCs should be greater than 1.\n")}
  if(maxdim>(nrow(D)-2)){stop("maxq should be less than the number of samples by at least two samples.\n")}
  if(missing(mindim)){mindim=2}
  if(missing(maxdim)){maxdim=3}
  is.wholenumber=function(x,tol=.Machine$double.eps^0.5) abs(x-round(x)) < tol
  if(!is.wholenumber(mindim)){stop("mindim should be a whole number.\n")}
  if(!is.wholenumber(maxdim)){stop("maxdim should be a whole number.\n")}
  if(mindim>maxdim){stop("mindim can not be greater than maxdim.\n")}
  if(SDselect<0){stop("SDselect can not be less than zero.\n")}
  if(maxdim>ncol(D)){stop("maxdim can not be greater than the number of variables.\n")}
  if(maxdim>10){cat("Warning! Model fitting may become very slow for maxdim > 10.\n\n")}
  if(any(is.na(Conf))){Conf=rep(1,nrow(D))}
  type=names(table(Conf))
  nt=length(type)
  conF=as.matrix(Conf)
  ####################################
  ## installing missing R packages
  installed<-installed.packages()[,1] 
  required<-c("mvtnorm","mclust","sva","stats","ggplot2","RColorBrewer",
              "rARPACK","RcppArmadillo","RcppEigen","Rcpp","parallel",
              "foreach","doParallel","doMC","compiler","devtools")
  toinstall<-required[!(required %in% installed)]
  if(length(toinstall) != 0){
    source("https://bioconductor.org/biocLite.R")
    biocLite(toinstall)
  }
  lapply(required, require, character.only = TRUE)
  requirfMM<-c("fMM","exploBATCHcolon","exploBATCHbreast")
  toinstallfMM<-requirfMM[!(requirfMM %in% installed)]
  if(length(toinstallfMM) != 0){
      install_github("syspremed/fMM")
      install_github("syspremed/exploBATCHbreast")
      install_github("syspremed/exploBATCHcolon")
  }
  lapply(requirfMM, require, character.only = TRUE)
  enableJIT(1)
  ## theme for ggplot
  theme=theme(strip.text.y = element_text(),#rotate strip text to horizontal
              axis.text = element_text(colour = "black", family="Arial", size=15),
              axis.title.x = element_text(colour = "black", family="Arial",face="bold", size=15),
              axis.title.y = element_text(colour = "black", family="Arial",face="bold", size=15),
              axis.text.x = element_text(colour = "black", family="Arial", size=15),#remove y axis labels
              legend.background = element_rect(fill = "white"),
              panel.grid.major = element_line(colour = "white"),
              legend.title=element_blank(),
              legend.text=element_text(colour="black",family="Arial",size=15),
              strip.text.x = element_text(colour = "black", family="Arial", size = 15),
              plot.title=element_text(size=15, family="Arial", face="bold",hjust=.5,margin=margin(b=2,unit="pt")),
              panel.border=element_rect(colour="black",size=2),
              legend.key=element_blank()) #format title
  ## colors
  if(any(batchCL==0)){batchCL=batchCL+1}
  grps<-names(table(batchCL))
  ngrps<-length(grps)
  if(ngrps<9){cGsub=brewer.pal(n=8,name="Dark2")[1:ngrps]
  }else{cGsub=c(1:ngrps)}
  ## creating folder for storing results and setting the working directory
  dir.create("exploBATCHresults", showWarnings = FALSE)
  dir.create("exploBATCHresults/pcaBeforeCorrection", showWarnings = FALSE)
  dir.create("exploBATCHresults/ppccaBeforeCorrection", showWarnings = FALSE)
  dir.create("exploBATCHresults/findBATCH", showWarnings = FALSE)
  dir.create("exploBATCHresults/Rworkspace", showWarnings = FALSE)
  dirp<-getwd()   
  setwd(paste0(dirp,"/exploBATCHresults"))
  ############# data preprocessing for PCA and PPCCA
  op <- options(warn = (-1)) # suppress warnings
  Yus=as.matrix(D)
  Yus=Yus[,apply(Yus,2,sd)>SDselect]
  cat(ncol(Yus), " features with SD > ",SDselect," were selected. \n\n")
  if(ncol(Yus)==0){stop("No features were selected, reduce value of SDselect.\n")}
  if((ncol(Yus)<20)&(SDselect>0)){stop("No need for reducing the number features using SDselect since the number features is already small (less than 20).\n")}
  Ys=scalingfn(Yus,type=scale)     # scaling
  muhat=colMeans(Ys)
  tYc=t(Ys)-muhat                  # centering
  Yc=t(tYc)
  N=nrow(Yc); p=ncol(Yc)
  S=fM(tYc,Yc)/N
  sumdiS=mdiag(S,p)
  ## EVD to initialize models
  tryCatch(if(p>1000){tp=1000;tmp=eigs_sym(S,k=tp,which="LM");tp=length(!is.na(tmp$val))
  }else{tp=p;tmp=eigen(S)},error=function(e) NULL)
  ############## pca before batch correction
  sink(file="undesired_output.txt")
  rs1<-prcomp(Yc,scale=FALSE,center=FALSE)  #pca
  pcaBeforeCorrection(rs1,grps,cGsub,batchCL,Conf,type,Yus,theme)
  save(rs1,file="Rworkspace/pcaBeforeBatchCorrection.RData")
  rm(D,rs1)
  gc()
  ####
  ############## ppcca before batch correction
  n.cores=detectCores()
  res1<-ppcca(Y=Ys,Yc,tYc,X=batchCL,S,sumdiS,muhat,tp,tmp,mindim,maxdim,n.cores=n.cores)      # run ppcca
  ppccaBeforeCorrection(res1,grps,cGsub,batchCL,Conf,type,Yus,theme)
  qopt=res1$q
  save(res1,file="Rworkspace/ppccaBeforeBatchCorrection.RData")
  sink()
  ####
  ############## detect batch and biological effect using ppcca
  cat("Detecting batch effect via findBATCH...\n\n")
  sink(file="undesired_output.txt")
  res2<-ppccajack(Ys,res1,n.cores)        # run ppcca via Jacknifing
  save(res2,file="Rworkspace/findBATCH.RData")
  rm(res1)
  gc()
  rerun=FALSE
  if(nt>1){ # biological effects to retain
    res1j<-ppcca(Y=Ys,Yc,tYc,X=conF,S,sumdiS,muhat,tp,tmp,qopt,qopt,n.cores=n.cores)
    rerun<-ppccajack(Ys,res1j,n.cores)
    rm(res1j)}
  LCI<-findBATCH(res2,qopt,nt,rerun,theme)
  sink()
  gc()
  ####
  if(sum(LCI)>0){ ## correct batch effect if its significant
    dir.create("correctComBat", showWarnings = FALSE)
    ############## correct batch effect using combat
    sink(file="undesired_output.txt")
    comres<-correctComBat(Yus,batchCL,grps,cGsub,Conf,theme,type)
    save(comres,file="Rworkspace/ComBat.RData")
    sink()
    ############# assess combat corrected data
    dir.create("assessComBat", showWarnings = FALSE)
    cat("Assessing for batch effect after ComBat...\n\n")
    ## preprocess combat corrected data
    sink(file="undesired_output.txt")
    Ycbs=scalingfn(t(comres),type=scale)     # scaling
    muhatcb=colMeans(Ycbs)
    tYcbc=t(Ycbs)-muhatcb                  # centering
    Ycbc=t(tYcbc)
    Scb=fM(tYcbc,Ycbc)/N
    sumdiScb=mdiag(Scb,p)
    ## EVD to initialize models
    tryCatch(if(p>1000){tpcb=1000;tmpcb=eigs_sym(Scb,k=tpcb,which="LM");tpcb=length(!is.na(tmpcb$val))
    }else{tpcb=p;tmpcb=eigen(Scb)},error=function(e) NULL)
    res1cbj<-ppcca(Y=Ycbs,Ycbc,tYcbc,X=batchCL,Scb,sumdiScb,muhatcb,tpcb,tmpcb,qopt,qopt,n.cores=n.cores)
    rerun1<-ppccajack(Ycbs,res1cbj,n.cores)       # assess batch effect on combat corrected data
    save(rerun1,file="Rworkspace/assessComBat.RData")
    rm(res1cbj)
    gc()
    rerun12=FALSE
    if(nt>1){
      res2cbj<-ppcca(Ycbs,Ycbc,tYcbc,conF,Scb,sumdiScb,muhatcb,tpcb,tmpcb,qopt,qopt,n.cores=n.cores)  # assess 
      rerun12<-ppccajack(Ycbs,res2cbj,n.cores)
      rm(res2cbj)}
    # assess biological effect on combat corrected data
    assessComBat(rerun1,qopt,nt,rerun12,theme)
    sink()
    rm(Ycbs,tYcbc,Ycbc,Scb,muhatcb,sumdiScb,rerun1,rerun12)
    gc()
    cat("findBATCH FINISHED assessing for batch effect after ComBat!! \n\n")
    ####
    if(method=="ppcca"){
        
      dir.create("correctBATCH",showWarnings=FALSE)
      ############## correct batch effect using ppcca
      sink(file="undesired_output.txt")
      mindim=maxdim=N-2
      if((N>250)&(p>400)){mindim=maxdim=200}
      if(maxdim>ncol(tmp$vec)){mindim=maxdim=ncol(tmp$vec)-2}
      designX=as.matrix(model.matrix(~factor(batchCL))[,-1])
            
      res11=NULL
      while((is.null(res11))&(maxdim>150)){
        tryCatch({
          res11 <- ppcca(Y = Ys, Yc, tYc, X = designX, S, sumdiS, muhat, tp, tmp, mindim, maxdim, n.cores = n.cores)
        },error=function(ex) {print(ex$message)})
        mindim<-maxdim<-(maxdim-1)
      }
      
      ppccaDat<-correctBatch(res11,designX,Yus,batchCL,ngrps,grps,Conf,nt,cGsub,comres,theme,type)
      save.image("Rworkspace/correctBATCH.RData")
      sink()
      rm(Yus,Ys,tYc,Yc,S,muhat,sumdiS,res11)
      gc()
      cat("correctBATCH FINISHED correcting for batch effect!! \n\n")
      ############# assess ppcca corrected data
      dir.create("assessCorrectBATCH", showWarnings=FALSE)
      cat("Assessing for batch effect after correctBATCH...\n\n")
      sink(file="undesired_output.txt")
      ## preprocess ppcca corrected data
      Ys=scalingfn(ppccaDat,type=scale)     # scaling 
      rm(ppccaDat)
      gc()
      muhat=colMeans(Ys)
      tYc=t(Ys)-muhat                  # centering
      Yc=t(tYc)
      S=fM(tYc,Yc)/N
      sumdiS=mdiag(S,p)
      ## EVD to initialize models
      tryCatch(if(p>1000){tp=1000;tmp=eigs_sym(S,k=tp,which="LM");tp=length(!is.na(tmp$val))
      }else{tp=p;tmp=eigen(S)},error=function(e) NULL)
      res2j<-ppcca(Y=Ys,Yc,tYc,X=batchCL,S,sumdiS,muhat,tp,tmp,qopt,qopt,n.cores=n.cores)
      rerun2<-ppccajack(Ys,res2j,n.cores)       # assess batch effect on ppcca corrected data
      save.image("Rworkspace/assessCorrectBATCH.RData")
      rm(res2j)
      gc()
      rerun22=FALSE
      if(nt>1){res22j<-ppcca(Y=Ys,Yc,tYc,X=conF,S,sumdiS,muhat,tp,tmp,qopt,qopt,n.cores=n.cores)
      rerun22<-ppccajack(Ys,res22j,n.cores)}
      assessCorrectBATCH(rerun2,qopt,nt,rerun22,theme)
      sink()
      cat("findBATCH FINISHED assessing for batch effect after correctBATCH!! \n\n")}}
  ####
  sessionInfo()
  setwd(dirp)
  cat("exploBATCH completed running!! Output in your currect working directory. \n\n")
  options(op) # reset the default value
}
