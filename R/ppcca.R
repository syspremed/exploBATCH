ppcca <-
function(Y,Yc,tYc,X,S,sumdiS,muhat,tp,tmp, minq=2, maxq=3,eps=0.8,n.cores){
  X=as.matrix(X)
  options(cores = n.cores)
  #cl=makeCluster(n.cores, type = "FORK")
  #registerDoParallel(cl, cores=n.cores)
  #if(n.cores==1) registerDoSEQ() else registerDoParallel(cores=n.cores) # multicore functionality
  if(n.cores==1) registerDoSEQ() else   registerDoMC() # multicore functionality
  #mcoptions<-list(preschedule=TRUE)
  cat(paste("\nUsing",n.cores,"cores for parallel processing. \n\n"))
  V=300
  repV=rep(0,V)
  ndim=dim(Y)
  N=ndim[1]
  p=ndim[2]
  Ip=diag(p)
  picon=(p/2)*log(2*pi)
  Vp=10
  C2p=p*3
  Np=N*p
  Nvp=Np+Vp+2
  nT=length(minq:maxq)
  if(nT>1){  ## run the ppca over different pcs
    IQ=lapply(1:maxq,diag)
    penalty=lapply(1:maxq,function(q) pen(q,p,N))
    resQ<-foreach(q=minq:maxq,.export=c('emQ','chol.inv','Aitkenf','mdiag','fM'),.packages=c('fMM','RcppArmadillo','RcppEigen','Rcpp')) %dopar%{
      emQ(q=q,Yc,tYc,S,tmp,tp,p,Ip,Np,Nvp,C2p,picon,IQ,penalty,eps,repV,V)}
    ## selection of optimal number of pcs
    BIC=unlist(resQ)
    ad=(minq-1); qopt=c(minq:maxq-ad)[BIC[(minq:maxq-ad)]==max(BIC[(minq:maxq-ad)])]
    q=qopt+ad; Iq=IQ[[q]]
    BIC=BIC[minq:maxq-ad]
    #par(mfrow=c(1,1),oma=c(.5,.5,1,1),mar=c(2.7,2.7,0,0),mgp=c(1.5,.5,0),xpd=FALSE)
    #plot(minq:maxq,BIC,type="b",xlab="q",ylab="BIC",col=4,col.lab="blue",lwd=2,cex=1.3,cex.lab=1.5)
    #abline(v=qopt+ad,col="red",lty=2,lwd=2)
    #box(lwd=2)
  }else{BIC= rep(NA, maxq)
    q=minq=maxq
    Iq=diag(q)}
  ## run the ppcca 
  L=ncol(X)
  Covars=standardize(X)
  Covars=rbind(rep(1,N),t(Covars))
  Xproj=crossprod(Covars,chol.inv(tcrossprod(Covars)))
  # initialize ppcca
  Sig=abs(sum(tmp$val[(q+1):tp])/(tp-q))
  W=tmp$vec[,1:q]
  iM=chol.inv(crossprod(W)+Sig*Iq)
  iMtW=tcrossprod(iM,W)
  u=fM(iMtW,tYc)
  iniAlpha<-lapply(1:q,function(i) iniA(i,u,Covars,L))
  Alpha=unlistfn(iniAlpha,L+1,byrow=TRUE)
  del=Alpha%*%Covars
  ll=lla=repV
  tol=eps+1; v=0
  while(tol>eps){
    v=v+1
    u=fM(iMtW,tYc)+iM%*%(Sig*del)
    tu=t(u)
    sumEuu=N*Sig*iM+u%*%tu
    Alpha=u%*%Xproj
    del=Alpha%*%Covars
    W=(tYc%*%tu)%*%chol.inv(sumEuu)
    ww=crossprod(W)
    MLESig=(N*sumdiS+sum(ww*sumEuu)-2*sum(fM(Yc,W)*tu))/Np
    Sig=(Np*MLESig+C2p)/Nvp
    SIG=Sig*Iq
    M=ww+SIG
    iM=chol.inv(M)
    iMtW=tcrossprod(iM,W)
    InvS=(Ip-W%*%iMtW)/Sig
    tYcc=tYc-W%*%del
    MhalDist=colSums(tYcc*fM(InvS,tYcc))
    logdetD=p*log(Sig)+log(det(M/Sig))
    ll[v]=sum(-picon-.5*logdetD-.5*MhalDist)
    converge=Aitkenf(ll,lla,v,eps)
    tol=converge[[1]]
    lla[v]=converge[[2]]
    if(v==V){tol=eps-1}}
  cat("PPCCA converged with the optimal number of pcs  = ", q, "\n\n")
  list(q=q,Sig=Sig,u=u,W=W,coefficients=Alpha,BIC=BIC,Covars=Covars,eps=eps,V=V,repV=repV,L=L,p=p,Ip=Ip,Iq=Iq,picon=picon,Vp=Vp,C2p=C2p,Nvp=Nvp)}
