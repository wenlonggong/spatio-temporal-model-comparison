data<-read.csv("data.csv")

# Load necessary packages
library(sp)
library(gstat)
library(ncdf4)

pred<-matrix(NA,ncol=5,nrow=4488)
Var<-matrix(NA,ncol=5,nrow=4488)
Metric<-matrix(NA,ncol=5,nrow=5)
LLP <- NULL
ULP <- NULL

for(i in 1:5){
print(i)
  data<-read.csv(paste("train",i,'.csv',sep = ''))[,-1]
  test<-read.csv(paste("test",i,'.csv',sep = ''))[,-1]
  
# Create spatiotemporal points
data$Date_Hour <- as.POSIXct(data$Date_Hour, format = "%Y-%m-%d %H:%M:%S")
df<-stConstruct(data,space=c('Longitude','Latitude'),time='Date_Hour',SpatialObj=SpatialPoints(data[,c('Longitude','Latitude')]))
STIDF<- as(df, "STIDF")
sepVgm <- vgmST("separable", 
                space=vgm(1,"Exp", 1, 0.5),
                time =vgm(1,"Exp", 1, 0.5),
                sill=150)


# Create variogram model
variogram.model <- variogram(logPM2.5 ~ lY_t_1+lY_1+lY_2+lY_3+lY_4+lY_5+lY_6+lY_7+lY_8+lY_9+lY_10, data = STIDF,tunit="hours")
variogram.model<-fit.StVariogram(variogram.model,sepVgm)

################### transform test data

test$Date_Hour <- as.POSIXct(test$Date_Hour, format = "%Y-%m-%d %H:%M:%S")
df_test<-stConstruct(test,space=c('Longitude','Latitude'),time='Date_Hour',SpatialObj=SpatialPoints(test[,c('Longitude','Latitude')]))
STIDF_test<- as(df_test, "STIDF")

# Universal Kriging
uk <- krigeST(logPM2.5 ~ lY_t_1+lY_1+lY_2+lY_3+lY_4+lY_5+lY_6+lY_7+lY_8+lY_9+lY_10, data = STIDF, newdata = STIDF_test, modelList  = variogram.model,computeVar = TRUE)


pred[,i] <- as.matrix(uk@data[,1])
Var[,i] <- as.matrix(uk@data[,2])
LL <-   pred [,i]+qnorm(0.025)*sqrt(Var[,i])*pred [,i]
UL <-    pred [,i]+qnorm(0.975)*sqrt(Var[,i])*pred [,i]

LLP <- c(LLP,LL)
ULP <- c(ULP,UL)


Metric[i,1] = sum(abs((pred[,i]-test[,14])))/sum(abs((pred[,i]+test[,14])))
Metric[i,2] = sqrt(mean((pred[,i] - test[,14])^2))
Metric[i,3] = mean(abs(pred[,i]  - test[,14]))
Metric[i,4] = cor(pred[,i],test[,14])
Metric[i,5] = mean(test[,14] < UL & test[,14] > LL)
print(i)
}
colnames(Metric)<-c("SAMPE","RMSE", "MAD","Cor","Cov")
write.csv(Metric,"Krig_Metric_5folds.csv")
colnames(pred)<-c("fold1","fold2","fold3","fold4","fold5")
write.csv(pred,"Krig_prediction.csv")

colnames(Var)<-c("fold1","fold2","fold3","fold4","fold5")
write.csv(Var,"Krig_Var.csv")

apply(Metric,2,mean)
spplot(spd)
spplot(uk["var1.pred"], sp.layout = list("sp.points", data))

############################
# Predicting Grid##
###########################
spat_pred_grid <- expand.grid(
  lon = seq(min(data$Longitude), max(data$Longitude), by = 0.1),
  lat = seq(min(data$Longitude), max(data$Longitude), by = 0.1)) %>%
  SpatialPoints(proj4string = CRS(proj4string(df)))
temp_pred_grid = as.POSIXct("2019-06-17")+3600*(00:5)
DE_pred <- STF(sp = spat_pred_grid, # spatial part
               time = temp_pred_grid) # temporal part