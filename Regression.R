rm(list = ls())
train<-read.csv("train.csv")
nrow=i
Regression<-lm(logPM2.5 ~ lY_t_1+lY_1+lY_2+lY_3+lY_4+lY_5+lY_6+lY_7+lY_8+lY_9+lY_10,data=train)
Prediction<-predict(Regression,newdata =test)
pred<-matrix(Prediction,ncol=1)
Reg_prediction<-matrix(NA,ncol=5,nrow=i)
colnames(Reg_prediction)<-c("fold1","fold2","fold3","fold4","fold5")
Metric<-matrix(NA,ncol = 5,nrow=5)
for (i in 1:5){
  train<-read.csv(paste("train",i,'.csv',sep = ''))[,-1]
  test<-read.csv(paste("test",i,'.csv',sep = ''))[,-1]
  
  Regression<-lm(logPM2.5 ~ lY_t_1+lY_1+lY_2+lY_3+lY_4+lY_5+lY_6+lY_7+lY_8+lY_9+lY_10+Krig,data=train)  
  Prediction<-predict(Regression,newdata =test,interval = "confidence")
  Reg_prediction[,i]<-matrix(Prediction[,1],nrow=i,ncol=1)
  
  Metric[i,1]<-sum(abs((Prediction[,1]-test[,14])))/sum(abs((Prediction[,1]+test[,14]))) 
  Metric[i,2] <- sqrt(mean((Prediction[,1] - test[,14])^2)) 
  Metric[i,3] = mean(abs(Prediction[,1] - test[,14]))  
  Metric[i,4] = cor(Prediction[,1],test[,14])            
  Metric[i,5]<- mean(test[,14]<Prediction[,3] & test[,14] >Prediction[,2]) 
}