ppccajack <-
function(Y,res2,n.cores){
  options(cores = n.cores)
  if(n.cores==1) registerDoSEQ() else   registerDoMC() # multicore functionality
  V=res2$V
  eps=res2$eps
  repV=res2$repV
  p=res2$p
  Ip=res2$Ip
  q=res2$q
  Iq=res2$Iq
  L=res2$L
  Vp=res2$Vp
  C2p=res2$C2p
  picon=res2$picon
  N=nrow(Y)
  nj=N-1
  Np=nj*p
  Nvp=Np+Vp+2
  Wopt=res2$W
  Sigopt=res2$Sig
  Covars=res2$Covars
  Alphaj=Alphaopt=res2$coefficients
  pnames=rownames(Covars)
  if(is.null(pnames)){
    pnames=c(1:(L+1))
    pnames[1]="Intercept"
    pnames[2:(L+1)]=paste("Covr",1:L,sep="-")}
  pnames[1]="Intercept"
  phi=1; Jack=ceiling(phi*N) #Jack=N
  idx=1:Jack
  cat("Starting jackniffing to capture parameter uncertainty...\n \n")
  resJ<-foreach(n=1:Jack,.export=c('jEM','chol.inv','Aitkenf','mdiag','fM','scalingfn','rowVar'),.packages=c('fMM','RcppArmadillo','RcppEigen','Rcpp')) %dopar%{
    jEM(n,Y,idx,nj,p,Covars,Alphaopt,Wopt,Sigopt,Iq,repV,eps,Np,Nvp,C2p,Ip,picon,V)}
  se_Alpha=sqrt((Jack-1)^2/Jack*(apply(simplify2array(resJ),c(1,2),var)))
  Ztab=qnorm(0.975)*se_Alpha
  UpCIAlp=Alphaopt+Ztab
  LoCIAlp=Alphaopt-Ztab
  colnames(LoCIAlp)=paste("LoCI-",pnames,sep="")
  colnames(UpCIAlp)=paste("UpCI-",pnames,sep="")
  tab1=cbind(LoCIAlp,UpCIAlp)
  for(i in 1:(L+1)){if(i==1){CI_Alp=cbind(tab1[,i],tab1[,(i+L+1)])
  colnames(CI_Alp)=colnames(tab1)[c(i,(i+L+1))]
  }else{tab2=cbind(tab1[,i], tab1[,(i+L+1)])
  colnames(tab2)=colnames(tab1)[c(i,(i+L+1))]
  CI_Alp=cbind(CI_Alp,tab2)}}
  list(coefficients=Alphaopt, coeffCI=CI_Alp)}
