rm(list = ls())
library(RcppZiggurat)
library(Rcpp)
library(Rfast)
library(GpGp)
library(tidyr)
library(graphics)
library(sp)
library(ggmap)
library( viridis)
PA<-read.csv("data.csv")[,-1]
PA_coor_unique<-unique(PA[,6:5])



loc<-PA[,c("Longitude","Latitude")]
locs<-as.matrix(loc)
response<-PA[,14]
X<-as.matrix( rep(1,nrow(locs)))
fit <- fit_model(response, locs, X,"exponential_isotropic",max_iter = 1000,start_parms=c(0.29545030, 0.03628063,0.01728957),convtol = 1e-05,reorder = TRUE)


distance.matrix <- dista(PA[,c('Longitude','Latitude')], PA[,c('Longitude','Latitude')],type = "euclidean")
rule=which(distance.matrix > 0.4815336, arr.ind=T )
#loc.rule<-PA_coor_unique[rule[,2],]


for(i in 1:748){
  PA_coor_unique[rule[1,2],]
}



### 10 Nearest Neighbor locations

dd<-matrix(NA)
data<-matrix(NA,nrow=1,ncol=41)
colnames(data)<-c("dd","Longitude", "Latitude","Y_1","T","Longitude", "Latitude","Y_2","T","Longitude","Latitude","Y_3","T","Longitude","Latitude","Y_4","T",
                  "Longitude", "Latitude","Y_5","T","Longitude", "Latitude","Y_6","T","Longitude","Latitude","Y_7","T","Longitude","Latitude","Y_8","T","Longitude","Latitude","Y_9","T","Longitude","Latitude","Y_10","T")
for(i in 1:935){
  dd<-matrix(NA)
  for(j in 1:10){
    target<-PA_coor_unique[order_dist_to_point(locs, locs[i,])[2:11],][j,1]
    try<-filter(Purpleair_2019,Longitude%in%target)[,c(6,5,12,13)]
    colnames(try)<-c("Longitude","Latitude",paste("Y",j,sep = '_'),"T")
    dd<-cbind(dd,try)
  }
  data<-rbind(data,dd)
}