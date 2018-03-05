setwd("/home/zhangl4/jinxx493/mpmri/")
library(polspline)
library(parallel)
library(spatstat)
library(ROCR)
library(MASS)
library(Matrix)
library(sparseMVN)
library(mvtnorm)
library(fields)
index=1:46
n.cores=100
load("newcasename.RData")


rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

matern.cor <- function(phi, nu, V){
  R=matrix(NA,nrow(V),nrow(V))
  for (ss in 1:nrow(V))
  {
    for (tt in 1:nrow(V))
    {
      R[ss,tt]=Matern(d=V[ss,tt],alpha=phi,nu=nu)
    }
  }
  return(R)
}


#### Using complete data to get region specific MRI distributions:
load("server_fillrecenter.RData")
fillin=fillrecenter
rm(fillrecenter)
for (i in 1:46)
{
  fillin[[i]][,"ADC"]=fillin[[i]][,"ADC"]/100
  fillin[[i]][,c("KTRANS","KEP","AUGC")]=log(fillin[[i]][,c("KTRANS","KEP","AUGC")])
}
########## Generate data for jags: 
data=cbind(as.data.frame(fillin[[1]][,c("ADC","KTRANS","KEP","AUGC","cancer","case","x","y","zone")]),1)
colnames(data)=c("ADC","KTRANS","KEP","AUGC","cancer","case","x","y","zone","subject")
for (k in c(index[-1]))
{
  temp=cbind(as.data.frame(fillin[[k]][,c("ADC","KTRANS","KEP","AUGC","cancer","case","x","y","zone")]),k)
  colnames(temp)=c("ADC","KTRANS","KEP","AUGC","cancer","case","x","y","zone","subject")
  data=rbind(data,temp)
}
dataNC=data[data[,"cancer"]==0,]
dataC=data[data[,"cancer"]==1,]
DATA=data
rm(data)
meancpz=colMeans(DATA[((DATA[,"cancer"]==1)&(DATA[,"zone"]=="PZ")),c("ADC","KTRANS","KEP","AUGC")])
meanccg=colMeans(DATA[((DATA[,"cancer"]==1)&(DATA[,"zone"]=="CG")),c("ADC","KTRANS","KEP","AUGC")])
meanncpz=colMeans(DATA[((DATA[,"cancer"]==0)&(DATA[,"zone"]=="PZ")),c("ADC","KTRANS","KEP","AUGC")])
meannccg=colMeans(DATA[((DATA[,"cancer"]==0)&(DATA[,"zone"]=="CG")),c("ADC","KTRANS","KEP","AUGC")])
means=list(meancpz,meanccg,meanncpz,meannccg)

covcpz=var(DATA[((DATA[,"cancer"]==1)&(DATA[,"zone"]=="PZ")),c("ADC","KTRANS","KEP","AUGC")])
covccg=var(DATA[((DATA[,"cancer"]==1)&(DATA[,"zone"]=="CG")),c("ADC","KTRANS","KEP","AUGC")])
covncpz=var(DATA[((DATA[,"cancer"]==0)&(DATA[,"zone"]=="PZ")),c("ADC","KTRANS","KEP","AUGC")])
covnccg=var(DATA[((DATA[,"cancer"]==0)&(DATA[,"zone"]=="CG")),c("ADC","KTRANS","KEP","AUGC")])
covs=list(covcpz,covccg,covncpz,covnccg)
####### average pcancer per region:
ppz=nrow(DATA[(DATA[,"cancer"]==1)&(DATA[,"zone"]=="PZ"),])/nrow(DATA[(DATA[,"zone"]=="PZ"),])
pcg=nrow(DATA[(DATA[,"cancer"]==1)&(DATA[,"zone"]=="CG"),])/nrow(DATA[(DATA[,"zone"]=="CG"),])
poverall=nrow(dataC)/nrow(DATA)

Means=list()

Means[[1]]=list()
Means[[1]][[1]]=means[[1]]
Means[[1]][[2]]=means[[2]]+c(-2,-0.5,-0.3,-0.3)
Means[[1]][[3]]=means[[3]]
Means[[1]][[4]]=means[[4]]+c(-2,-0.5,-0.3,-0.3)

Means[[2]]=means

Means[[3]]=list()
Means[[3]][[1]]=means[[1]]
Means[[3]][[2]]=means[[1]]
Means[[3]][[3]]=means[[3]]
Means[[3]][[4]]=means[[3]]

Means[[4]]=list()
Means[[4]][[1]]=means[[1]]
Means[[4]][[2]]=means[[2]]+c(-2,-0.5,-0.3,-0.3)
Means[[4]][[3]]=means[[3]]+c(2,-0.5,-0.5,-0.3)
Means[[4]][[4]]=means[[4]]+c(2,-0.5,-0.5,-0.3)+c(-2,-0.5,-0.3,-0.3)

Means[[5]]=list()
Means[[5]][[1]]=means[[1]]
Means[[5]][[2]]=means[[2]]
Means[[5]][[3]]=means[[3]]+c(2,-0.5,-0.5,-0.3)
Means[[5]][[4]]=means[[4]]+c(2,-0.5,-0.5,-0.3)

Means[[6]]=list()
Means[[6]][[1]]=means[[1]]
Means[[6]][[2]]=means[[1]]
Means[[6]][[3]]=means[[3]]+c(2,-0.5,-0.5,-0.3)
Means[[6]][[4]]=means[[3]]+c(2,-0.5,-0.5,-0.3)

#the1=c(0.4,0.8)  #sigmasq
#the2=c(3,7) #phi
#the3=c(1,2,3,4,5,6)
#setting=expand.grid(x=the1,y=the2,z=the3)
 
the1=c(0.4,0.8)  #sigmasq
the2=c(5,10) #phi
the3=c(1,2,3)
setting1=expand.grid(x=the1,y=the2,z=the3)
the1=c(0.2,0.8)  #sigmasq
the2=c(10,15) #phi
the3=c(1,2,3)
setting2=expand.grid(x=the1,y=the2,z=the3)
setting=rbind(setting1,setting2)
 
index=1:44
set=1:1
  sigmasq=setting[set,1]
  #phi=setting[set,2]
  phi=runif(9,3,12)
  nu=setting[set,3]
  means=Means[[nu]]
  auc_mbase1=tp_mbase1=numeric()
  auc_mregion1=tp_mregion1=numeric()
  auc_mcoord1=tp_mcoord1=numeric()
  auc_msmooth1=tp_msmooth1=numeric()
  
  auc_mbase=tp_mbase=numeric()
  auc_mregion=tp_mregion=numeric()
  auc_mcoord=tp_mcoord=numeric()
  auc_msmooth=tp_msmooth=numeric()
  
results_mbase=results_mregion=results_mcoord=results_msmooth=list()
for (s in 1:10)
{
  results_mbase[[s]]=list()
  results_mregion[[s]]=list()
  results_mcoord[[s]]=list()
  results_msmooth[[s]]=list()
############ Simulate Data:
##### get index:
  select=sample(1:46,size=44,replace=TRUE)
  load("server_fillrecenter_missing.RData")
  fillin=fillrecenter
  rm(fillrecenter)
  for (i in 1:46)
  {
    fillin[[i]][,"ADC"]=fillin[[i]][,"ADC"]/100
    fillin[[i]][,c("KTRANS","KEP","AUGC")]=log(fillin[[i]][,c("KTRANS","KEP","AUGC")])
  }
  
  ########## Generate data for jags: 
  data=cbind(as.data.frame(fillin[[1]][,c("ADC","KTRANS","KEP","AUGC","cancer","case","x","y","zone")]),1)
  colnames(data)=c("ADC","KTRANS","KEP","AUGC","cancer","case","x","y","zone","subject")
  for (k in c(index[-1]))
  {
    temp=cbind(as.data.frame(fillin[[k]][,c("ADC","KTRANS","KEP","AUGC","cancer","case","x","y","zone")]),k)
    colnames(temp)=c("ADC","KTRANS","KEP","AUGC","cancer","case","x","y","zone","subject")
    data=rbind(data,temp)
  }
  dataNC=data[data[,"cancer"]==0,]
  dataC=data[data[,"cancer"]==1,]
  DATA=data
  rm(data)
  
  simulate_fillnew=function(i)
  {
    temp=fillin[[i]]
    zone=1*(temp[,"zone"]=="PZ")+0*(temp[,"zone"]=="CG")
    #rannorm = sapply(1:nrow(temp),FUN=function(x){rnorm(1, 0, sqrt(tausq))})
    coords=temp[,c("x","y")]
    locations=cbind(cords,1)
    for (rows in 1:nrow(coords))
    {
      locations[rows,3]=1*((coords[rows,1]<=(-1/3))&(coords[rows,2]<=(-1/3))) 
      + 2*((coords[rows,1]<=(-1/3))&(coords[rows,2]<=(1/3))&(coords[rows,2]>(-1/3)))
      + 3*((coords[rows,1]<=(-1/3))&(coords[rows,2]>(1/3)))
      + 4*((coords[rows,1]<=(1/3))&(coords[rows,1]>(-1/3))&(coords[rows,2]<(-1/3)))
      + 5*((coords[rows,1]<=(1/3))&(coords[rows,1]>(-1/3))&(coords[rows,2]<=(1/3))&(coords[rows,2]>(-1/3)))
      + 6*((coords[rows,1]<=(1/3))&(coords[rows,1]>(-1/3))&(coords[rows,2]>(1/3)))
      + 7*((coords[rows,1]>(1/3))&(coords[rows,2]<=(-1/3)))
      + 8*((coords[rows,1]>(1/3))&(coords[rows,2]<=(1/3))&(coords[rows,2]>(-1/3)))
      + 9*((coords[rows,1]>(1/3))&(coords[rows,2]>(1/3)))
    }
    D <- as.matrix(dist(coords))
    R=D
    for (rows in 1:nrow(coords))
    {
      for (cols in 1:nrow(coords))
      {
        R[rows,cols]=1/2*exp(-phi[locations[rows,3]]*D[rows,cols]) + 1/2*exp(-phi[locations[cols,3]]*D[rows,cols])
      }
    }
    #qs <- mvrnorm(n=1,mu=rep(0,nrow(temp)),Sigma=sigmasq*R)+sapply(1:nrow(temp),FUN=function(x){qnorm(p=zone[x]*ppz+(1-zone[x])*pcg,mean=0,sd=1)})
    qs <- rmvn(1, rep(0,nrow(temp)), sigmasq*R)+sapply(1:nrow(temp),FUN=function(x){qnorm(p=zone[x]*ppz+(1-zone[x])*pcg,mean=0,sd=1)})
    while (sum(qs>0)==0) {
      qs <- rmvn(1, rep(0,nrow(temp)), sigmasq*R)+sapply(1:nrow(temp),FUN=function(x){qnorm(p=zone[x]*ppz+(1-zone[x])*pcg,mean=0,sd=1)})
    }
    ps <- sapply(1:nrow(temp),FUN=function(x){pnorm(q=qs[x],mean=0,sd=1)})
    #cancer=sapply(1:nrow(temp),FUN=function(x){rbinom(n=1,size=1,prob=pc[x])})
    cancer=1*(ps>poverall)+0*(ps<=poverall)
    temp[,"cancer"]=cancer
    temp=cbind(temp,ps,qs)
    for (j in 1:nrow(temp))
    {
      indicator=1*((cancer[j]==1)&(zone[j]==1))+2*((cancer[j]==1)&(zone[j]==0))+3*((cancer[j]==0)&(zone[j]==1))+4*((cancer[j]==0)&(zone[j]==0))
      temp[j,c("ADC","KTRANS","KEP","AUGC")]=mvrnorm(n=1,mu=means[[indicator]],Sigma=covs[[indicator]])
    }
    return(temp)
  }
  
  fillnew=mclapply(1:44,FUN=function(x){simulate_fillnew(x)},mc.cores=44)
  
fillin=fillnew
index=1:44

for(i in 1:44)
{
  temp=fillin[[i]]
  temp[,"area"]=1*(temp[,"cancer"]==1)+2*((temp[,"cancer"]==0)&(temp[,"zone"]=="PZ"))+3*((temp[,"cancer"]==0)&(temp[,"zone"]=="CG"))
  temp[,"parea"]=1*(temp[,"cancer"]==1)+2*((temp[,"cancer"]==0)&(temp[,"Z"]==2))+3*((temp[,"cancer"]==0)&(temp[,"Z"]==1))
  temp[,"region"]=1*((temp[,"cancer"]==1)&(temp[,"zone"]=="PZ"))+2*((temp[,"cancer"]==1)&(temp[,"zone"]=="CG"))+3*((temp[,"cancer"]==0)&(temp[,"zone"]=="PZ"))+4*((temp[,"cancer"]==0)&(temp[,"zone"]=="CG"))
  temp[,"pregion"]=1*((temp[,"cancer"]==1)&(temp[,"Z"]==2))+2*((temp[,"cancer"]==1)&(temp[,"Z"]==1))+3*((temp[,"cancer"]==0)&(temp[,"Z"]==2))+4*((temp[,"cancer"]==0)&(temp[,"Z"]==1))
  fillin[[i]]=temp
}
fillrecenter=fillin
fillnew=fillin
########## Generate data by slice: ########## 
################ Unsmoothed noimpute:
##### Redefine!
##### Redefine!
fillrecenter[[1]][,"subject"]=rep(1,nrow(fillrecenter[[1]]))
data=fillrecenter[[1]]
for (k in 2:44)
{
  #A=as.data.frame(cbind(fillrecenter[[k]],k))
  #colnames(A)=c(names(fillrecenter[[1]]),"subject")
  fillrecenter[[k]][,"subject"]=rep(k,nrow(fillrecenter[[k]]))
  data=rbind(data,fillrecenter[[k]])
}
#data[,c("KTRANS","KEP","AUGC")]=log(data[,c("KTRANS","KEP","AUGC")])
dataC=data[data[,"cancer"]==1,]
dataNC=data[data[,"cancer"]==0,]
DATA=data;DATAC=dataC;DATANC=dataNC
rm(data)

rdata=list()
ldata=list()
ydata=list()
adata=list()
padata=list()
region=list()
pregion=list()
zone=list()
TZ=list()
Z=list()
ZP=list()

for (K in 1:44)
{
  rdata[[K]]=cbind(fillrecenter[[K]][,c("ADC","KTRANS","KEP","AUGC")])
  #rdata[[K]][,c("KTRANS","KEP","AUGC")]=log(rdata[[K]][,c("KTRANS","KEP","AUGC")])
  ldata[[K]]=as.matrix(cbind(fillrecenter[[K]][,c("|x|","y","zone")]))
  colnames(ldata[[K]])=c("loc.x","loc.y","loc")
  ydata[[K]]=fillrecenter[[K]][,c("cancer")]
  adata[[K]]=fillrecenter[[K]][,c("area")]
  padata[[K]]=fillrecenter[[K]][,c("parea")]
  region[[K]]=fillrecenter[[K]][,c("region")]
  pregion[[K]]=fillrecenter[[K]][,c("pregion")]
  zone[[K]]=fillrecenter[[K]][,c("zone")]
  Z[[K]]=fillrecenter[[K]][,c("Z")]
  TZ[[K]]=fillrecenter[[K]][,c("TZ")]
  ZP[[K]]=fillrecenter[[K]][,c("ZP")]
}

k=1
n.new=1
llist =matrix(0:1,nrow=2)

alpha=0.5
m=4 ################################## # of parameters!!
delta=m-1 #
p=4
npat=1


######### CV to tune smoothing bandwith:
if (s==1)
{
B=seq(0.04,0.5,by=0.02)
source("SBDA_spline_CancerX.R") 
PRED=function(data.pred, data.train, alpha){
  k=1
  dec.mat=array(list(NULL), c(npat, 3))
  logpden=apply(llist, 1, FUN=function(x){LogPredDen(data.train, data.pred[k,,drop=F], x, delta,omega=omega.hat)})
  pden=exp(logpden-mean(logpden))
  
  ##### P
  PHAT=data.pred[[1,"paxis"]]
  
  priprob=apply(llist, 1, FUN=function(x){k=sum(x==1); PHAT^k*(1-PHAT)^(n.new-k)})
  pden=priprob*pden
  #####
  pden=pden/sum(pden)
  
  risk=sapply(1:nrow(llist) , FUN=function(loop){
    d=llist[loop,]
    n.fp = apply(llist, 1,FUN=function(y){sum(d==1 & y==0)})
    n.fn = apply(llist, 1,FUN=function(y){sum(d==0 & y==1)})
    risk=sum((pden)*((1-alpha)*n.fn+alpha*n.fp))
  })
  
  dec.mr=as.numeric(llist[which.min(risk),])
  dec.mat[[k,1]]=data.pred[[k,4]]
  dec.mat[[k,2]]=dec.mr
  dec.mat[[k,3]]=pden
  return(dec.mat)
}
whole=DATA[,-1]
wholePZ=DATA[DATA$zone=="PZ",-1]
wholeCG=DATA[DATA$zone=="CG",-1]

bpre=btru=Spre=list()
ss=1
for (j in 1:44)
{
  bpre[[ss]]=numeric()
  Spre[[ss]]=numeric()
  btru[[ss]]=numeric()
  
  ##### trainng:
  trainC=DATAC[DATAC[,"subject"]!=j,-1]
  trainNC=DATANC[DATANC[,"subject"]!=j,-1]
  train=rbind(trainC,trainNC)
  trainPZ=DATA[((DATA[,"subject"]!=j)&(DATA[,"zone"]=="PZ")),-1]
  trainCG=DATA[((DATA[,"subject"]!=j)&(DATA[,"zone"]=="CG")),-1]
  
  fitNC=polymars(trainNC[,c("ADC","KTRANS","KEP","AUGC")], trainNC[,c("x","y")], classify = F,maxsize=6,knots=9)
  fitC=polymars(trainC[,c("ADC","KTRANS","KEP","AUGC")], trainC[,c("x","y")], classify = F,maxsize=6,knots=9)
  
  X=list()
  for (qq in 1:44)
  {
    X[[qq]]=list()
    temp=whole[whole[,"subject"]==qq,]
    X[[qq]][[1]]=design.polymars(fitNC,temp[,c("x","y")])
    X[[qq]][[2]]=design.polymars(fitC,temp[,c("x","y")])
  }
  dim_beta=c(ncol(X[[1]][[1]]),ncol(X[[1]][[2]]))
  
  ######## FIT P 
  fitp=polymars(train[,"cancer"], train[,c("|x|","y")], classify = T,maxsize=4)
  W=list()
  for (qq in 1:44)
  {
    temp=whole[whole[,"subject"]==qq,]
    W[[qq]]=predict(fitp,temp[,c("|x|","y")],classify=F)[,1]
  }
  
  ############## Training Data:
  mydata.train=data.frame()
  mydata.train[1:43,1]=index[-j]
  mydata.train$par=rdata[index[-j]]
  mydata.train$axis=X[index[-j]]
  mydata.train$cancer=ydata[index[-j]]
  mydata.train$paxis=W[index[-j]]
  mydata.train$area=adata[index[-j]]
  mydata.train$region=region[index[-j]]
  mydata.train$zone=zone[index[-j]]
  mydata.train$predzone=Z[index[-j]]
  mydata.train$predzonep=ZP[index[-j]]
  
  #phat=sum(unlist(mydata.train[,5])==1)/length(unlist(mydata.train[,5]))
  dec.mat=array(list(NULL), c(npat, 3))
  
  ParmFit=DataFit(mydata.train)
  Nreg=omega.hat=NULL
  for(z in 1:2)
  {omega.hat[z]=list(ParmFit$om.hat[[z]]*delta/(length(ParmFit$r[[z]])-1))}
  
  ###### Spline X for beta:
  
  simu=function(P)
  {
    mydata.test=data.frame()
    mydata.test[1,1]=j
    mydata.test$par=list(rdata[[j]][P,])
    mydata.test$axis=list(list(X[[j]][[1]][P,],X[[j]][[2]][P,]))
    mydata.test$cancer=list(ydata[[j]][P])
    mydata.test$paxis=list(W[[j]][P])
    mydata.test$area=list(adata[[j]][P])
    mydata.test$region=list(region[[j]][P])
    mydata.test$zone=list(zone[[j]][P])
    mydata.test$predzone=list(Z[[j]][P])
    mydata.test$predzonep=list(ZP[[j]][P])
    a=PRED(data.pred=mydata.test,data.train=mydata.train, alpha)
    return(a)
  }
  
  results=mclapply(1:nrow(rdata[[j]]),simu,mc.cores = n.cores)
  
  for (i in 1:length(results))
  {
    bpre[[ss]][i]=as.numeric(results[[i]][1,3][[1]][2]/sum(results[[i]][1,3][[1]]))
    btru[[ss]][i]=as.numeric(results[[i]][1,1][[1]])
  }
  #bpred <- prediction(bpre[[ss]], btru[[ss]])
  #bperf<- performance(bpred,"tpr","fpr")
  #bau <- performance(bpred,"auc")
  #auc[j] <- unlist(slot(bau, "y.values"))
  #tp[j]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
  ss=ss+1
  #print(ss)
}
#tuning:
auc=tp=list()
ff=1
for (bb in B)
{
  auc[[ff]]=tp[[ff]]=numeric()
for (j in 1:44)
{
  ###### msmooth:
  Xcoord=fillin[[j]][,"x"];Ycoord=fillin[[j]][,"y"];
  XX <- ppp(Xcoord, Ycoord, c(-1,1),c(-1,1),marks=bpre[[j]])
  #b <- bw.smoothppp(X)  # least-squares CV to select a smoothing bandwidth for spatial smoothing of marks.
  Spre=as.numeric(Smooth(XX, sigma=bb,at="points",edge=TRUE, diggle=FALSE))
  #### Msmooth:
  bpred <- prediction(Spre, btru[[j]])
  bperf<- performance(bpred,"tpr","fpr")
  bau <- performance(bpred,"auc")
  auc[[ff]][j] <- unlist(slot(bau, "y.values"))
  tp[[ff]][j]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
}
  ff=ff+1
  #print(ff)
}  

AUC=sapply(1:length(B),FUN=function(x){mean(auc[[x]])})
b=B[which.max(AUC)]
}


####### Mbase:
source("SBDA_F.R")

PRED=function(data.pred, data.train, alpha,omega,ParmFit){
  k=1
  dec.mat=array(list(NULL), c(npat, 3))
  logpden=apply(llist, 1, FUN=function(x){LogPredDen(data.train, data.pred[k,,drop=F], x, delta,omega,ParmFit)})
  pden=exp(logpden-mean(logpden))
  
  priprob=apply(llist, 1, FUN=function(x){k=sum(x==1); phat^k*(1-phat)^(n.new-k)})
  pden=priprob*pden
  
  risk=sapply(1:nrow(llist) , FUN=function(loop){
    d=llist[loop,]
    n.fp = apply(llist, 1,FUN=function(y){sum(d==1 & y==0)})
    n.fn = apply(llist, 1,FUN=function(y){sum(d==0 & y==1)})
    risk=sum((pden)*((1-alpha)*n.fn+alpha*n.fp))
  })
  
  dec.mr=as.numeric(llist[which.min(risk),])
  dec.mat[[k,1]]=data.pred[[k,4]]
  dec.mat[[k,2]]=dec.mr
  dec.mat[[k,3]]=pden
  return(dec.mat)
}

############################ Mbase:
#### Training:
mydata.train=data.frame()
mydata.train[1:34,1]=1:34
mydata.train$par=rdata[1:34]
mydata.train$axis=ldata[1:34]
mydata.train$cancer=ydata[1:34]

phat=sum(unlist(mydata.train[,4])==1)/length(unlist(mydata.train[,4]))
dec.mat=array(list(NULL), c(npat, 3))

ParmFit=DataFit(mydata.train)
Nreg=omega.hat=NULL
for(z in 1:2)
{omega.hat[z]=list(ParmFit$sigma.hat[[z]]*delta/(length(ParmFit$r[[z]])-1))}

bpre=btru=list()
ss=1
auc=tp=numeric()
for (j in 35:44)
{
  simu=function(P)
  {
    mydata.test=data.frame()
    mydata.test[1,1]=j
    mydata.test$par=list(rdata[[j]][P,])
    mydata.test$axis=list(ldata[[j]][P,])
    mydata.test$cancer=list(ydata[[j]][P])
    a=PRED(data.pred=mydata.test,data.train=mydata.train, alpha,omega=omega.hat,ParmFit)
    return(a)
  }
  results=mclapply(1:nrow(rdata[[j]]),simu,mc.cores=n.cores)
  results_mbase[[s]][[ss]]=results
  bpre[[ss]]=numeric()
  btru[[ss]]=numeric()
  for (i in 1:length(results))
  {
    bpre[[ss]][i]=results[[i]][1,3][[1]][2]/sum(results[[i]][1,3][[1]])
    btru[[ss]][i]=as.numeric(results[[i]][1,1][[1]])
  }
  bpred <- prediction(bpre[[ss]], btru[[ss]])
  bperf<- performance(bpred,"tpr","fpr")
  bau <- performance(bpred,"auc")
  auc[ss] <- unlist(slot(bau, "y.values"))
  tp[ss]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
  ss=ss+1
}
auc_mbase1[s]=mean(auc)
tp_mbase1[s]=mean(tp)
bpre=unlist(bpre)
btru=unlist(btru)
  bpred <- prediction(bpre, btru)
  bperf<- performance(bpred,"tpr","fpr")
  bau <- performance(bpred,"auc")
  auc_mbase[s] <- unlist(slot(bau, "y.values"))
  tp_mbase[s]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
  print(s)
  print(auc_mbase1[s])

#################### Mregion:
source("4group.R")
PRED=function(data.pred, data.train, alpha){
  k=1
  dec.mat=array(list(NULL), c(npat, 3))
  logpden=apply(llist, 1, FUN=function(x){LogPredDen(data.train, data.pred[k,,drop=F], x, delta,omega=omega.hat)})
  pden=exp(logpden-mean(logpden))
  
  ##### P
  PHAT=data.pred[[1,3]]
  priprob=apply(llist, 1, FUN=function(x){k=sum(x==1); PHAT^k*(1-PHAT)^(n.new-k)})
  pden=priprob*pden
  
  risk=sapply(1:nrow(llist) , FUN=function(loop){
    d=llist[loop,]
    n.fp = apply(llist, 1,FUN=function(y){sum(d==1 & y==0)})
    n.fn = apply(llist, 1,FUN=function(y){sum(d==0 & y==1)})
    risk=sum((pden)*((1-alpha)*n.fn+alpha*n.fp))
  })
  
  dec.mr=as.numeric(llist[which.min(risk),])
  dec.mat[[k,1]]=data.pred[[k,4]]
  dec.mat[[k,2]]=dec.mr
  dec.mat[[k,3]]=pden
  return(dec.mat)
}


n.new=1
llist =matrix(0:1,nrow=2)
spauc=numeric()

whole=DATA[,-1]
wholePZ=DATA[DATA$zone=="PZ",]
wholeCG=DATA[DATA$zone=="CG",]

###### Spline X for beta:
trainC=DATAC[(DATAC[,"subject"]%in%c(1:34)),]
trainNC=DATANC[(DATANC[,"subject"]%in%c(1:34)),]
trainPZ=DATA[((DATA[,"subject"]%in%c(1:34))&(DATA$zone=="PZ")),]
trainCG=DATA[((DATA[,"subject"]%in%c(1:34))&(DATA$zone=="CG")),]
train=DATA[(DATA[,"subject"]%in%c(1:34)),]

test=DATA[(DATA[,"subject"]%in%c(35:44)),]
testNC=DATANC[(DATANC[,"subject"]%in%c(35:44)),]

### PZ:2, CG:1
train[,"ZONE"]=ifelse(train[,"zone"]=="PZ",2,1)
fitp=polymars(train[,"cancer"], as.factor(train[,"ZONE"]), classify = T,maxsize=4)
W=list()
for (qq in 1:44)
{
  temp=whole[whole[,"subject"]==qq,]
  p0=predict(fitp,as.factor(temp[,"Z"]),classify=F)[,1]
  OW=ifelse(temp[,"Z"]==2,1,2)
  p1=predict(fitp,as.factor(OW),classify=F)[,1]
  pPZ=ifelse(temp[,"Z"]==2,p0,p1)
  pCG=ifelse(temp[,"Z"]==1,p0,p1)
  W[[qq]]=pPZ*ZP[[qq]]+pCG*(1-ZP[[qq]])
}

############## Training Data:
mydata.train=data.frame()
mydata.train[1:34,1]=1:34
mydata.train$par=rdata[1:34]
mydata.train$paxis=W[1:34] # predicted p
mydata.train$cancer=ydata[1:34]
mydata.train$area=adata[1:34]
mydata.train$region=region[1:34]
mydata.train$predzone=Z[1:34]  #### dat[[i,7]]: zone for test voxel! in logpredden
mydata.train$predzonep=ZP[1:34]

phat=sum(unlist(mydata.train[,4])==1)/length(unlist(mydata.train[,4]))
dec.mat=array(list(NULL), c(npat, 3))

ParmFit=DataFit(mydata.train)
Nreg=omega.hat=NULL
for(z in 1:4)
{omega.hat[z]=list(ParmFit$sigma.hat[[z]]*delta/(length(ParmFit$r[[z]])-1))}

bpre=btru=list()
auc=tp=numeric()
ss=1
for (j in 35:44)
{
  bpre[[ss]]=numeric()
  btru[[ss]]=numeric()
    simu=function(P)
    {
      mydata.test=data.frame()
      mydata.test[1,1]=j
      mydata.test$par=list(rdata[[j]][P,])
      mydata.test$paxis=list(W[[j]][P])
      mydata.test$cancer=list(ydata[[j]][P])
      mydata.test$area=list(adata[[j]][P])
      mydata.test$region=list(region[[j]][P])
      mydata.test$predzone=list(Z[[j]][P])
      mydata.test$predzonep=list(ZP[[j]][P])
      a=PRED(data.pred=mydata.test,data.train=mydata.train, alpha)
      return(a)
    }
    results=mclapply(1:nrow(rdata[[j]]),simu,mc.cores = n.cores)
    results_mregion[[s]][[ss]]=results
    for (i in 1:length(results))
    {
      bpre[[ss]][i]=as.numeric(results[[i]][1,3][[1]][2]/sum(results[[i]][1,3][[1]]))
      btru[[ss]][i]=as.numeric(results[[i]][1,1][[1]])
    }
    bpred <- prediction(bpre[[ss]], btru[[ss]])
    bperf<- performance(bpred,"tpr","fpr")
    bau <- performance(bpred,"auc")
    auc[ss] <- unlist(slot(bau, "y.values"))
    tp[ss]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
    ss=ss+1
}
auc_mregion1[s]=mean(auc)
tp_mregion1[s]=mean(tp)
    bpre=unlist(bpre)
    btru=unlist(btru)
    bpred <- prediction(bpre, btru)
    bperf<- performance(bpred,"tpr","fpr")
    bau <- performance(bpred,"auc")
    auc_mregion[s] <- unlist(slot(bau, "y.values"))
    tp_mregion[s]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
    print(auc_mregion1[s])

######################## Mcoord & Msmooth:
source("SBDA_spline_CancerX.R") 
PRED=function(data.pred, data.train, alpha){
  k=1
  dec.mat=array(list(NULL), c(npat, 3))
  logpden=apply(llist, 1, FUN=function(x){LogPredDen(data.train, data.pred[k,,drop=F], x, delta,omega=omega.hat)})
  pden=exp(logpden-mean(logpden))
  
  ##### P
  PHAT=data.pred[[1,"paxis"]]
  
  priprob=apply(llist, 1, FUN=function(x){k=sum(x==1); PHAT^k*(1-PHAT)^(n.new-k)})
  pden=priprob*pden
  #####
  pden=pden/sum(pden)
  
  risk=sapply(1:nrow(llist) , FUN=function(loop){
    d=llist[loop,]
    n.fp = apply(llist, 1,FUN=function(y){sum(d==1 & y==0)})
    n.fn = apply(llist, 1,FUN=function(y){sum(d==0 & y==1)})
    risk=sum((pden)*((1-alpha)*n.fn+alpha*n.fp))
  })
  
  dec.mr=as.numeric(llist[which.min(risk),])
  dec.mat[[k,1]]=data.pred[[k,4]]
  dec.mat[[k,2]]=dec.mr
  dec.mat[[k,3]]=pden
  return(dec.mat)
}
whole=DATA[,-1]
wholePZ=DATA[DATA$zone=="PZ",-1]
wholeCG=DATA[DATA$zone=="CG",-1]

trainC=DATAC[DATAC[,"subject"]%in%c(1:34),-1]
trainNC=DATANC[DATANC[,"subject"]%in%c(1:34),-1]
train=rbind(trainC,trainNC)
trainPZ=DATA[((DATA[,"subject"]%in%c(1:34))&(DATA[,"zone"]=="PZ")),-1]
trainCG=DATA[((DATA[,"subject"]%in%c(1:34))&(DATA[,"zone"]=="CG")),-1]

fitNC=polymars(trainNC[,c("ADC","KTRANS","KEP","AUGC")], trainNC[,c("x","y")], classify = F,maxsize=6,knots=9)
fitC=polymars(trainC[,c("ADC","KTRANS","KEP","AUGC")], trainC[,c("x","y")], classify = F,maxsize=6,knots=9)

X=list()
for (qq in 1:44)
{
  X[[qq]]=list()
  temp=whole[whole[,"subject"]==qq,]
  X[[qq]][[1]]=design.polymars(fitNC,temp[,c("x","y")])
  X[[qq]][[2]]=design.polymars(fitC,temp[,c("x","y")])
}
dim_beta=c(ncol(X[[1]][[1]]),ncol(X[[1]][[2]]))

######## FIT P 
fitp=polymars(train[,"cancer"], train[,c("|x|","y")], classify = T,maxsize=4)
W=list()
for (qq in 1:44)
{
  temp=whole[whole[,"subject"]==qq,]
  W[[qq]]=predict(fitp,temp[,c("|x|","y")],classify=F)[,1]
}

############## Training Data:
mydata.train=data.frame()
mydata.train[1:34,1]=1:34
mydata.train$par=rdata[1:34]
mydata.train$axis=X[1:34]
mydata.train$cancer=ydata[1:34]
mydata.train$paxis=W[1:34]
mydata.train$area=adata[1:34]
mydata.train$region=region[1:34]
mydata.train$zone=zone[1:34]
mydata.train$predzone=Z[1:34]
mydata.train$predzonep=ZP[1:34]

#phat=sum(unlist(mydata.train[,5])==1)/length(unlist(mydata.train[,5]))
dec.mat=array(list(NULL), c(npat, 3))

ParmFit=DataFit(mydata.train)
Nreg=omega.hat=NULL
for(z in 1:2)
{omega.hat[z]=list(ParmFit$om.hat[[z]]*delta/(length(ParmFit$r[[z]])-1))}

bpre=btru=Spre=list()
auc=tp=numeric()
sauc=stp=numeric()
ss=1
for (j in 35:44)
{
  
  bpre[[ss]]=numeric()
  Spre[[ss]]=numeric()
  btru[[ss]]=numeric()
    ###### Spline X for beta:
    
    simu=function(P)
    {
      mydata.test=data.frame()
      mydata.test[1,1]=j
      mydata.test$par=list(rdata[[j]][P,])
      mydata.test$axis=list(list(X[[j]][[1]][P,],X[[j]][[2]][P,]))
      mydata.test$cancer=list(ydata[[j]][P])
      mydata.test$paxis=list(W[[j]][P])
      mydata.test$area=list(adata[[j]][P])
      mydata.test$region=list(region[[j]][P])
      mydata.test$zone=list(zone[[j]][P])
      mydata.test$predzone=list(Z[[j]][P])
      mydata.test$predzonep=list(ZP[[j]][P])
      a=PRED(data.pred=mydata.test,data.train=mydata.train, alpha)
      return(a)
    }
    
    results=mclapply(1:nrow(rdata[[j]]),simu,mc.cores = n.cores)
    results_mcoord[[s]][[ss]]=results
    for (i in 1:length(results))
    {
      bpre[[ss]][i]=as.numeric(results[[i]][1,3][[1]][2]/sum(results[[i]][1,3][[1]]))
      btru[[ss]][i]=as.numeric(results[[i]][1,1][[1]])
    }
    bpred <- prediction(bpre[[ss]], btru[[ss]])
    bperf<- performance(bpred,"tpr","fpr")
    bau <- performance(bpred,"auc")
    auc[ss] <- unlist(slot(bau, "y.values"))
    tp[ss]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
    
    ###### msmooth:
    Xcoord=fillin[[j]][,"x"];Ycoord=fillin[[j]][,"y"];
    XX <- ppp(Xcoord, Ycoord, c(-1,1),c(-1,1),marks=bpre[[ss]])
     #b <- bw.smoothppp(X)  # least-squares CV to select a smoothing bandwidth for spatial smoothing of marks.
    Spre[[ss]]=as.numeric(Smooth(XX, sigma=b,at="points",edge=TRUE, diggle=FALSE))
    bpred <- prediction(Spre[[ss]], btru[[ss]])
    bperf<- performance(bpred,"tpr","fpr")
    bau <- performance(bpred,"auc")
    sauc[ss] <- unlist(slot(bau, "y.values"))
    stp[ss]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
    ss=ss+1
}
#### Mcoord:
    auc_mcoord1[s]=mean(auc)
    tp_mcoord1[s]=mean(tp)
    bpre=unlist(bpre)
    btru=unlist(btru)
    bpred <- prediction(bpre, btru)
    bperf<- performance(bpred,"tpr","fpr")
    bau <- performance(bpred,"auc")
    auc_mcoord[s] <- unlist(slot(bau, "y.values"))
    tp_mcoord[s]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
    print(auc_mcoord1[s])

#### Msmooth:
    auc_msmooth1[s]=mean(sauc)
    tp_msmooth1[s]=mean(stp)
    Spre=unlist(Spre)
    bpred <- prediction(Spre, btru)
    bperf<- performance(bpred,"tpr","fpr")
    bau <- performance(bpred,"auc")
    auc_msmooth[s] <- unlist(slot(bau, "y.values"))
    tp_msmooth[s]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
    print(auc_msmooth1[s])
    
    save(auc_mbase,tp_mbase,auc_mregion,tp_mregion,auc_mbase1,tp_mbase1,auc_mregion1,tp_mregion1,
         auc_mcoord,tp_mcoord,auc_msmooth,tp_msmooth,auc_mcoord1,tp_mcoord1,auc_msmooth1,tp_msmooth1,
         b, results_mbase,results_mregion,results_mcoord,
         file=paste("simulation_paper1/simu_nonstationary_sigmasq_",sigmasq,"_nu",nu,".RData",sep="")) 
}
