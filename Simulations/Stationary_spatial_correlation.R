############ When predicting X (design matrix), the 
source("~/Google Drive/mpMRI/1004/SBDA_spline_CancerX_predmissing.R") ####### 2, without zone

PRED=function(data.pred, data.train, alpha){
  k=1
  dec.mat=array(list(NULL), c(npat, 3))
  #if (sum(is.na(data.pred[1,2][[1]]))==4){pden=c(1,1)}
  #if (sum(is.na(data.pred[1,2][[1]]))<4)
  #{
  logpden=apply(llist, 1, FUN=function(x){LogPredDen(data.train, data.pred[1,,drop=F], x, delta,omega=omega.hat)})
  pden=exp(logpden-mean(logpden))
  #}
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

missing_best_s2=matrix(NA,34,length(maxsize))
result_missing_best_s2=dec.mat=array(list(NULL), c(11, 34))

whole=DATA[,-1]
wholePZ=DATA[DATA$zone=="PZ",-1]
wholeCG=DATA[DATA$zone=="CG",-1]
whole2=DATAT[,-1]
wholePZ2=DATAT[DATAT$zone=="PZ",-1]
wholeCG2=DATAT[DATAT$zone=="CG",-1]


alpha=0.5
m=4 # of parameters!!
delta=m-1 #
p=4
npat=1
n.new=1
llist =matrix(0:1,nrow=2)
tpr=numeric()
spc=numeric()
ppv=numeric()
npv=numeric()
true=list()
pred=list()
spauc=numeric()
maxsize=5:15

M=1 ##### Change this!!
tp_msmooth_missing=numeric()
INDEX=1:length(unique(filename[[2]]))
INDEX.name=unique(filename[[2]])
for (k in 1:34)
{
  u=k
  for (j in u)
  {
  ###### Spline X for beta:
    trainC=DATAC[DATAC[,"subject"]!=j,-1]
    trainNC=DATANC[DATANC[,"subject"]!=j,-1]
    train=rbind(trainC,trainNC)
    testNC=DATANC[DATANC[,"subject"]==j,-1]
    trainPZ=DATA[((DATA[,"subject"]!=j)&(DATA[,"zone"]=="PZ")),-1]
    trainCG=DATA[((DATA[,"subject"]!=j)&(DATA[,"zone"]=="CG")),-1]
    
  fitNC=polymars(trainNC[,c("ADC","KTRANS","KEP","AUGC")], trainNC[,c("x","y")], classify = F,maxsize=6,knots=10)
  fitC=polymars(trainC[,c("ADC","KTRANS","KEP","AUGC")], trainC[,c("x","y")], classify = F,maxsize=6,knots=10)
  
  X=list()
  for (qq in index[-j])
  {
    X[[qq]]=list()
    temp=whole[whole[,"subject"]==qq,]  ### Also include missing when predicting.
    X[[qq]][[1]]=design.polymars(fitNC,temp[,c("x","y")])
    X[[qq]][[2]]=design.polymars(fitC,temp[,c("x","y")])
  }
  for (qq in index[j])
  {
    X[[qq]]=list()
    temp=whole2[whole2[,"subject"]==qq,]  ### Also include missing when predicting.
    X[[qq]][[1]]=design.polymars(fitNC,temp[,c("x","y")])
    X[[qq]][[2]]=design.polymars(fitC,temp[,c("x","y")])
  }
  dim_beta=c(ncol(X[[1]][[1]]),ncol(X[[1]][[2]]))
  
  ######## FIT P 
  fitp=polymars(train[,"cancer"], train[,c("|x|","y")], classify = T,maxsize=4)
  W=list()
  for (qq in index[-j])
  {
    temp=whole[whole[,"subject"]==qq,]
    W[[qq]]=predict(fitp,temp[,c("|x|","y")],classify=F)[,1]
  }
  for (qq in index[j])
  {
    temp=whole2[whole2[,"subject"]==qq,]
    W[[qq]]=predict(fitp,temp[,c("|x|","y")],classify=F)[,1]
  }
  
  ############## Training Data:
  mydata.train=data.frame()
  mydata.train[1:(34-length(u)),1]=index[-u]
  mydata.train$par=rdata[-u]
  mydata.train$axis=X[-u]
  mydata.train$cancer=ydata[-u]
  mydata.train$paxis=W[-u]
  mydata.train$area=adata[-u]
  mydata.train$region=region[-u]
  mydata.train$zone=zone[-u]
  mydata.train$predzone=Z[-u]
  mydata.train$predzonep=ZP[-u]
  
  #phat=sum(unlist(mydata.train[,5])==1)/length(unlist(mydata.train[,5]))
  dec.mat=array(list(NULL), c(npat, 3))
  
  ParmFit=DataFit(mydata.train)
  Nreg=omega.hat=NULL
  for(z in 1:2)
  {omega.hat[z]=list(ParmFit$om.hat[[z]]*delta/(length(ParmFit$r[[z]])-1))}
  
  simu=function(P)
  {
    mydata.test=data.frame()
    mydata.test[1,1]=index[j]
    mydata.test$par=list(rdata2[[j]][P,])
    mydata.test$axis=list(list(X[[j]][[1]][P,],X[[j]][[2]][P,]))
    mydata.test$cancer=list(ydata2[[j]][P])
    mydata.test$paxis=list(W[[j]][P])
    mydata.test$area=list(adata2[[j]][P])
    mydata.test$region=list(region2[[j]][P])
    mydata.test$zone=list(zone2[[j]][P])
    mydata.test$predzone=list(Z2[[j]][P])
    mydata.test$predzonep=list(ZP2[[j]][P])
    a=PRED(data.pred=mydata.test,data.train=mydata.train, alpha)
    #a=PRED_missing(data.pred=mydata.test,data.train=mydata.train, alpha)
    return(a)
  }
  
  results=mclapply(1:nrow(rdata2[[j]]),simu,mc.cores = 4)
  result_missing_best_s2[[M,j]]=results
  
  bpre=numeric()
  btru=numeric()
  for (i in 1:length(results))
  {
    bpre[i]=as.numeric(results[[i]][1,3][[1]][2]/sum(results[[i]][1,3][[1]]))
    btru[i]=as.numeric(results[[i]][1,1][[1]])
  }
  
  Xcoord=fillrecenter[[j]][,"x"];Ycoord=fillrecenter[[j]][,"y"];
  X <- ppp(Xcoord, Ycoord, c(-1,1),c(-1,1),marks=bpre)
  b=0.2
  #b <- bw.smoothppp(X)  # least-squares CV to select a smoothing bandwidth for spatial smoothing of marks.
  Spre=as.numeric(Smooth(X, sigma=b,at="points",edge=TRUE, diggle=F))
  
  bpred <- prediction(Spre, btru)
  bperf<- performance(bpred,"tpr","fpr")
  bau <- performance(bpred,"auc")
  
  missing_best_s2[j,M] <- unlist(slot(bau, "y.values"))
  tp_msmooth_missing[j]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
  print(j)
  print(missing_best_s2[j,M])
  }
}

save(result_missing_best_s2,missing_best_s2,file="~/Google Drive/mpMRI/1004/missing_best_s2_34.RData") 

#colnames(missing_best_s2)=c("xy12-16_|x|y2","xy15-13_|x|y2","xy15-13_|x|y4","xy14-14_|x|y4",rep(0,16))
#### Best: 14-14-4 0.825 11.13.

