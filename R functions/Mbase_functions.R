loghconst=function(delta, Omega){
  delta/2*as.numeric(determinant(Omega, logarithm = TRUE)$modulus)-0.5*nrow(Omega)*delta*log(2)- logGamma_p(0.5*delta, p=nrow(Omega))
}

logGamma_p=function(a, p){
  p*(p-1)/4*log(pi)+sum(lgamma(a+(1-1:p)/2))
}

###################################################################################
# DataFit: Estimate the adjusted mean y.tilde and covariance S.tile in (7) and (8)#
###################################################################################

each=function(i){
  keep=dat[[i,4]]==(z-1)
  n.i=sum(keep)
  ###################
  Y.i=as.matrix((dat[[i,2]])[keep,]) 
  if(n.i==1){
    Y.i=matrix(as.numeric(Y.i), nrow=1)
  }
  if(n.i==0){
    y.tilde.i=rep(0, p) 
    S.tilde.i=array(0, dim=c(p,p))
    phidet.i=NULL
    r.i=NULL
  }
  if (n.i>=1)
  {
    r.i=rep(1,n.i)
    y.tilde.i=apply(r.i*Y.i, 2, sum)
    S.tilde.i= t(Y.i)%*%Y.i
    phidet.i=1 
  }
  y.tilde=y.tilde+y.tilde.i 
  S.tilde=S.tilde+S.tilde.i
  r.all=c(r.all, r.i)	
  phidet=c(phidet, phidet.i)
  return()
}

DataFit=function(dat){
  #################################################################################
  ##  return y.tilde and S.tilde
  ## correlation function implied by variable struc  specified at the beginning
  #################################################################################
  
  Npat=nrow(dat)
  p=ncol(dat[[1,2]])
  mu.hat=vector(mode="list", 2)
  sigma.hat=vector(mode="list", 2)
  phidet.hat=vector(mode="list", 2)
  r=vector(mode="list", 2)
  
  for(z in 1:2){
    r.all=phidet=NULL
    y.tilde=rep(0, p)
    S.tilde=array(0,dim=c(p,p))
   
    for (i in 1:Npat)
    {
      keep=dat[[i,4]]==(z-1)
      n.i=sum(keep)
      ###################
      if(n.i>0){Y.i=as.matrix((dat[[i,2]])[keep,])}
      if(n.i==0){
        y.tilde.i=rep(0, p) 
        S.tilde.i=array(0, dim=c(p,p))
        phidet.i=NULL
        r.i=NULL
      }
      if (n.i>=1)
      {
        y.tilde.i=colSums(Y.i)
        S.tilde.i= t(Y.i)%*%Y.i
        phidet.i=1 
      }
      y.tilde=y.tilde+y.tilde.i 
      S.tilde=S.tilde+S.tilde.i
      r.all=c(r.all, rep(1,n.i))	
      phidet=c(phidet, phidet.i)
    }
    y.tilde=y.tilde/sum(r.all)
    S.tilde=S.tilde-y.tilde%*%t(y.tilde)*sum(r.all)
    mu.hat[z]=list(y.tilde)
    sigma.hat[z]=list(S.tilde)
    phidet.hat[z]=list(phidet)
    r[z]=list(r.all)
  }
  output=list(mu.hat=mu.hat, sigma.hat=sigma.hat, phidet=phidet.hat,r=r)
  return(output)
}


LogPredDen=function(data.train, data.pred, d, delta, omega, Parmfit){
  stopifnot(nrow(data.pred[[2]])==length(d))
  
  data.pred[[1,4]]=as.numeric(d)
  data <- rbind(data.train, data.pred)
  #ParmFit <- DataFit(data)
  if (data.pred[[1,4]]==1) {z1=1;z2=0}
  if (data.pred[[1,4]]==0) {z1=0;z2=1}
  r.1=ParmFit$r[[z1+1]] 
  r.2=ParmFit$r[[z2+1]] 
  ytilde.z=(ParmFit$mu.hat[[z1+1]]*length(r.1)+data.pred[[1,2]][data.pred[[1,4]]==z1,])/(length(r.1)+1)
  Stilde.z=ParmFit$sigma.hat[[z1+1]]+(as.numeric(data.pred[[1,2]][data.pred[[1,4]]==z1,]))%*%t(as.numeric(data.pred[[1,2]][data.pred[[1,4]]==z1,]))+(length(r.1))*ParmFit$mu.hat[[z1+1]]%*%t(ParmFit$mu.hat[[z1+1]])-(length(r.1)+1)*as.numeric(ytilde.z)%*%t(as.numeric(ytilde.z))
  output=-p/2*(log(sum(r.2)))-loghconst(delta+length(r.2)-1,  ParmFit$sigma.hat[[z2+1]]+omega[[z2+1]])-p/2*(log(sum(r.1)+1))-loghconst(delta+length(r.1)+1-1,  Stilde.z+omega[[z1+1]])
  return(output)
}
##################################################################
# MaxLLrho: maximum marginal likelihood given the correlation rho#
##################################################################

MarginLLrho=function(data, rho, omega){
  ### marginal likelihood p(Y|rho, omega)
  output=0
  Npat=nrow(data)
  ParmFit=DataFit(data,rho)
  for(z in 0:1)
  {
    Stilde.z=ParmFit$sigma.hat[[z+1]]
    omega.z=omega[[z+1]]
    psidet.z=ParmFit$phidet[[z+1]]
    r.z=ParmFit$r[[z+1]]
    output=output-p/2*sum(log(psidet.z))-p/2*log(sum(r.z))-loghconst(delta+length(r.z)-1,  Stilde.z+omega.z)
  }
  return(output)
}


MaxLLrho=function(data, rho){
  ### marginal likelihood p(Y|rho) with maximized omega.hat
  ParmFit=DataFit(data,rho)
  out=0
  for(z in 0:1){
    Stilde.z=ParmFit$sigma.hat[[z+1]]
    psidet.z=ParmFit$phidet[[z+1]]
    r.z=ParmFit$r[[z+1]]
    omega.z=Stilde.z*delta/(length(r.z)-1) ## (18)
    det1=as.numeric(determinant(Stilde.z+omega.z, logarithm = TRUE)$modulus)
    det2=as.numeric(determinant(omega.z, logarithm = TRUE)$modulus)
    out=out+(-p/2*sum(log(psidet.z))-p/2*log(sum(r.z))+loghconst(delta, omega.z)-loghconst(delta+length(r.z)-1, omega.z+Stilde.z))
  }
  return(out)
}

listToMat=function(data){
  dataMat=NULL
  for(i in 1:nrow(data)){
    dataMat=rbind(dataMat, data[[i,2]])
  }
  return(dataMat)
}
