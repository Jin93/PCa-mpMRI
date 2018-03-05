setwd("/home/reillyc/jinxx493/mpmri")
library(ROCR)
library(polspline)
library(parallel)
library(spatstat)
n.cores=34
index=1:46
load("server_fillrecenter.RData")
FILL=fillrecenter

Data=as.data.frame(cbind(fillrecenter[[1]],1))
colnames(Data)=c(names(fillrecenter[[1]]),"subject")
for (I in c(2:46))
{
  A=as.data.frame(cbind(fillrecenter[[I]],I))
  colnames(A)=c(names(fillrecenter[[1]]),"subject")
  Data=rbind(Data,A)
}
Data[,c("KTRANS","KEP","AUGC")]=log(Data[,c("KTRANS","KEP","AUGC")])
DataC=Data[Data[,"cancer"]==1,]
DataNC=Data[Data[,"cancer"]==0,]
DATA=Data;DATAC=DataC;DATANC=DataNC
Data0=DATA
load("newcasename.RData")
INDEX=1:length(unique(filename[[2]]))
INDEX.name=unique(filename[[2]])

index=1:34

source("CancerX_best_parallel_34.R")
load("newindex_34.RData")

k=1
n.new=1
llist =matrix(0:1,nrow=2)

#auc_mbase.boot=list()
#tp_mbase.boot=list()
#meanauc_mbase.boot=meantp_mbase.boot=numeric()
alpha=0.5
m=4 ################################## # of parameters!!
delta=m-1 #
p=4
npat=1

results_msmooth_bootstrap=mclapply(1:2,FUN=function(x){boot.best(x)},mc.cores=25)
save(results_msmooth_bootstrap,file="results_msmooth_bootstrap.RData")

