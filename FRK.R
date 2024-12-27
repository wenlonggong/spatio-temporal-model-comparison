library(anytime)
# getting the data
#data<-extract_Purpleair_12
data<-read.csv("data.csv")
data$Date_Hour<-ymd_hms(data$Date_Hour)
data<-subset(data,data$t1>=90 & data$t1<92) #just two days 164 obs


BAUs1 <- auto_BAUs(STplane(),type = "grid",data = STIDF,cellsize = c(0.01, 0.01, 1),convex=-0.4,tunit = "hours")    
BAUs2 <- auto_BAUs(STplane(),type = "grid",data = new,convex=-0.6,tunit = "hours")    




  
STIDF$std <- 2 # fix the measurement error variance
G_spatial <- auto_basis(manifold = plane(),      # fns on plane
                        data = as(STIDF, "Spatial"), # project
                        nres = 2,                     # 2 res.
                        type = "bisquare",            # bisquare.
                        regular = 0)                  # irregular

## ------------------------------------------------------------------------
t_grid <- matrix(seq(1, 48, length = 48))

## ------------------------------------------------------------------------
G_temporal <- local_basis(manifold = real_line(),  # fns on R1
                          type = "bisquare",       # bisquare
                          loc = t_grid,            # centroids
                          scale = rep(2, 48))      # aperture par.

## ------------------------------------------------------------------------
G <- TensorP(G_spatial, G_temporal)      # take the tensor product

## ------------------------------------------------------------------------
g1 <- show_basis(G_spatial) + xlab("lon (deg)") + ylab("lat (deg)") + coord_fixed()
g2 <- show_basis(G_temporal) + xlab("t (hours)") +ylab(expression(varphi(t)))
g1;g2;
STIDF$std <- sqrt(0.049)

## ------------------------------------------------------------------------
f <- logPM2.5 ~ 1

## ------------------------------------------------------------------------
S <- FRK(f = f,               # formula
         data = list(STIDF), # (list of) data
         basis = G,           # basis functions
         BAUs = BAUs1,         # BAUs
         n_EM = 3,            # max. no. of EM iterations
         tol = 0.01)          # tol. on change in log-likelihood

## ------------------------------------------------------------------------
locs_pred<-matrix(c(-118.2269,34.06887), 48, 2,byrow=T)
locs_pred<-matrix(c(-118.1586,34.06659), 48, 2,byrow=T)
colnames(locs_pred)<-c('Longitude','Latitude')

D<-as.data.frame(cbind(locs_pred,two_days[,c(8,7,29)]))
two_days<-subset(l_12,l_12$t>=90 & l_12$t<92)
two_days$Longitude<-locs_pred[,1]
two_days$Latitude<-locs_pred[,2]
newdata<- stConstruct(two_days,space=c('Longitude','Latitude'),crs = CRS("+proj=longlat +ellps=WGS84"),time='DATE_TIME')
coordinates(two_days) = ~ Longitude+Latitude
proj4string(two_days) = "+proj=longlat +datum=WGS84"
anytime(as.factor("2013-06-01 08:07:00"))
STIDF_pred<- as(pred_grid, "STFDF")
STIDF_pred<- as(newdata, "STIDF")
grid_BAUs <- predict(S,pred_grid)
  
stplot(grid_BAUs,main = "Predictions (degrees Fahrenheit)",layout = c(3, 2))
##
sp = SpatialPoints(locs_pred)
time = two_days$datetime

stidf = as(STIDF(sp, time, two_days), "STIDF")
stidf[1:2,]
all.equal(stidf, stidf[stidf,])

####
grid_BAUs@sp  <- SpatialPoints(grid_BAUs@sp)   # convert to Spatial points object
gridded(grid_BAUs@sp) <- TRUE                          # and assert that it is gridded
library(RColorBrewer)
colour_regions <- (brewer.pal(11, "Spectral") %>%     # construct spectral brewer palette
                     colorRampPalette())(16) %>%    # at 16 levels
  rev()                         # reverse direction

grid_BAUs$mu_trunc <- pmin(pmax(grid_BAUs$mu,68),105)
stplot(grid_BAUs[,c(4,9,14,19,24,29),"mu"],             # plot the FRK predictor
       main="Predictions (degrees Fahrenheit)",          # title
       layout=c(3,2),                                     # trellis layout
       col.regions=colour_regions,                        # color scale
       #xlim=c(-100,-80),ylim=c(32,46),                    # axes limits
       aspect=1)  
stplot(grid_BAUs[,c(4,9,14,19,24,29),"sd"],             # plot the FRK predictor
       main="Predictions (degrees Fahrenheit)",          # title
       layout=c(3,2),                                     # trellis layout
       col.regions=colour_regions,                        # color scale
       #xlim=c(-100,-80),ylim=c(32,46),                    # axes limits
       aspect=1)  




pred_grid1 <- expand.grid(
  lon =-118.2269,
  lat = 34.06659,hour=1:48)
pred_grid2 <- expand.grid(
  lon =seq(-118.24,-118.3,length=1),
  lat=seq(34.1,34.1,length=1),hour=1:48)
pred_grid<-rbind(pred_grid1,pred_grid2)
ll<-subset(l_12,l_12$t>=90 & l_12$t<92)
pred_grid$Time<-rep(ymd_hms(ll$DATE_TIME),1)
pred_grid<-pred_grid[,-3]
new<-stConstruct(pred_grid,space=c('lon','lat'),crs = CRS("+proj=longlat +ellps=WGS84"),time='Time')
new<- as(new, "STIDF")
## ------------------------------------------------------------------------
S_pred <- eval_basis(basis = G,                    # basis functs
                     s = pred_grid[,c("lon","lat","hour")] %>% # pred locs
                       as.matrix())%>%            # conv. to matrix
  as.matrix()                                 # as matrix
colnames(S_pred) <- paste0("B", 1:ncol(S_pred))  # assign  names
pred_grid <- cbind(pred_grid,S_pred)             # attach to grid
# options(error=recover)
# reach_full_in <- reachability(krack_full, 'in')
## ------------------------------------------------------------------------
pred_EPA12 <- predict(S,BAUs2)
























data("NOAA_df_1990")
Tmax <- subset(NOAA_df_1990, month %in% 7 & year == 1993)
Tmax <- within(Tmax, {time = as.Date(paste(year,month,day,sep="-"))})
STObj <- stConstruct(x = Tmax, space = c("lon","lat"), time = "time", interval = TRUE)

## BAUs: spatial BAUs are 1x1 pixels, temporal BAUs are 1 day intervals
BAUs <- auto_BAUs(manifold = STplane(), 
                  cellsize = c(1, 1, 1),    
                  data=STObj, tunit = "days")
BAUs$fs <- 1 # scalar fine-scale variance matrix, implicit in previous examples

## Basis functions
G <- auto_basis(manifold = STplane(), data = STObj, nres = 2, tunit = "days")

## Run FRK
STObj$std <- 2 # fix the measurement error variance
S <- FRK(f = z ~ 1 + lat, data = list(STObj), 
         basis = G, BAUs = BAUs, est_error = FALSE, method = "TMB")
pred <- predict(S, percentiles = NULL)

