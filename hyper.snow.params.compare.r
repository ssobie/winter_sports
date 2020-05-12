##Snow-17 accumulation and ablation model.

source('/storage/home/ssobie/code/repos/winter_sports/site.coordinates.r',chdir=TRUE)

library(ncdf4)
library(sp)
##library(gstat)
library(fields)

##Based on Anderson (2006) and Mark Raleigh's matlab code.
##Primary Citations:
## 1.  Anderson, E. A. (1973), National Weather Service River Forecast System
## Snow   Accumulation   and   Ablation   Model,   NOAA   Tech.   Memo.   NWS
## HYDro-17, 217 pp., U.S. Dep. of Commer., Silver Spring, Md.
## 2.  Anderson, E. A. (1976), A point energy and mass balance model of a snow
## cover, NOAA Tech. Rep. 19, 150 pp., U.S. Dep. of Commer., Silver Spring, Md.


##-------------------------------------------------------------------------------
##Create netcdfs for kriged parameters

create_kriged_netcdf <- function(var.name,lon,lat,param,model,write.dir) {

  lon.atts <- list(standard_name="longitude",
                   long_name = "longitude",
                   units = "degrees_east",
                   axis = "X")

  lat.atts <- list(standard_name="latitude",
                   long_name = "latitude",
                   units = "degrees_north",
                   axis = "Y")

  scaling.atts <- list(standard_name = "Scaling Factor",
                       long_name = "Scaling Factor",
                       missing_value = -9999.0,
                       units = "")
  slope.atts <- list(standard_name = "Slope Factor",
                       long_name = "Slope Factor",
                       missing_value = -9999.0,
                       units = "")
  freq.atts <- list(standard_name = "Half Frequency Point",
                       long_name = "Half Frequency Point",
                       missing_value = -9999.0,
                       units = "")
  var.atts <- switch(var.name,
                     scale=scaling.atts,
                     slope=slope.atts,
                     freq=freq.atts)

  write.clim.name <- paste0(write.dir,var.name,'_hyper_snow_calibrated_parameter_',model,'_prism_TPS.nc')

  n.lon <- length(lon)
  n.lat <- length(lat)

  x.geog <- ncdim_def('lon', 'degrees_east', lon)
  y.geog <- ncdim_def('lat', 'degrees_north', lat)
  t.geog <- ncdim_def('time', 'days since 2000-01-01', 1,
                      unlim=FALSE, calendar='365_day')
  var.geog <- ncvar_def(var.name, units='m', dim=list(x.geog, y.geog, t.geog),
                        missval=-9999.0)
  file.nc <- nc_create(write.clim.name, var.geog)
  ncatt_put(file.nc,'time','standard_name','Time')
  ncatt_put(file.nc,'time','long_name','Time')

  print('Lon names')
  lon.names <- names(lon.atts)
  for (j in 1:length(lon.atts))
    ncatt_put(file.nc,varid='lon',attname=lon.names[j],attval=lon.atts[[j]])
  print('Lat names')
  lat.names <- names(lat.atts)
  for (j in 1:length(lat.atts))
    ncatt_put(file.nc,varid='lat',attname=lat.names[j],attval=lat.atts[[j]])
  print('Var names')
  var.names <- names(var.atts)
  for (j in 1:length(var.atts))
    ncatt_put(file.nc,varid=var.name,attname=var.names[j],attval=var.atts[[j]])

  global.atts <- list(institution="Pacific Climate Impacts Consortium (PCIC), Victoria, BC, www.pacificclimate.org",
                   contact="Pacific Climate Impacts Consortium",
                   Conventions="CF-1.4",
                   institute_id ="PCIC")
  global.names <- names(global.atts)
  for (g in 1:length(global.atts))
    ncatt_put(file.nc,varid=0,attname=global.names[g],attval=global.atts[[g]])

  ncvar_put(file.nc,'lon',lon)
  ncvar_put(file.nc,'lat',lat)
  ncvar_put(file.nc,var.name,param)

  nc_close(file.nc)
  return(write.clim.name)
}



get_add_coordinates <- function(site) {

   coordinates <- list(add_curd_mountain=c(-121.907555,50.107700,1000),
                      add_scuzzy_mountain=c(-121.632935,49.856295,1000),
                      add_russell_mountain=c(-122.049435,50.435671,1000),
                      add_hemlock_valley=c(-121.929363,49.383619,1000),
                      add_snowcap_peak=c(-122.655731,49.834085,1000),   
                      add_overseer_mountain=c(-123.440696,50.493472,1000),
                      add_richmond=c(-123.097766,49.165711,1000),
                      add_surrey=c(-122.833068,49.134959,1000))

  rv <- coordinates[[site]]
  return(rv)
}

get_parameters <- function(site) {
   site.params <- list(callaghan=c(46.25,0.875,2.8,1.0209),
                    grouse_mountain=c(51.5,0.81,3.3,1.0209),
                    orchid_lake=c(51.8,0.7,4.1,1.0209),
                    palisade_lake=c(51.6,0.825,3.15,1.0209),
                    dog_mountain=c(50.5,0.725,3.03,1.0209),
                    stave_lake=c(50.6,0.825,3.4,1.0209),
                    nahatlatch=c(45.45,0.8,3.08,1.0209),
                    wahleach=c(40.05,0.2,2.75,1.0209),
                    klesilkwa=c(40.1,0.175,2.67,1.0209),
                    brookmere=c(38,0.05,2.7,1.0209),
                    shovelnose_mountain=c(38.15,0.075,2.7,1.0209),
                    lightning_lake=c(40.15,0.05,2.75,1.0209),
                    hamilton_hill=c(38.45,0.1,2.85,1.0209),
                    upper_squamish=c(51.7,0.575,3.7,1.0209),
                    tenquille_lake=c(43.45,0.5,2.25,1.0209),
                    chilliwack_river=c(51.8,0.8,4.0,1.0209),
                    spuzzum_creek=c(51.8,0.875,4.1,1.0209),
                    blackwall_peak_pillow=c(45,0.5,2.8,1.0209),
                    wahleach_lake=c(51.75,0.825,4.1,1.0209),
                    add_curd_mountain=c(53.0,0.87,4,1.0209),
                    add_scuzzy_mountain=c(51.6,0.87,4.0,1.0209),
                    add_russell_mountain=c(51.6,0.87,4.0,1.0209),
                    add_hemlock_valley=c(44.0,0.6,3.0,1.0209),
                    add_snowcap_peak=c(51.6,0.87,4.0,1.0209),
                    add_overseer_mountain=c(53,0.87,4,1.0209),
                    add_richmond=c(50.6,0.8,4.0,1.0209),
                    add_surrey=c(50.6,0.8,4.0,1.0209))
   rv <- site.params[[site]]
   return(rv)
}

get_cal_parameters <- function(site,model) {

   cal.file <- paste0('/storage/data/projects/rci/data/winter_sports/calibration_parameters/',site,'_',model,'_more_limited.csv')
   cal.data <- read.csv(cal.file,header=T,as.is=T)[[1]]
   return(cal.data)
}

##---------------------------------------------------------------------
model <- 'PNWNAmet'

sites <- c('blackwall_peak_pillow','chilliwack_river','spuzzum_creek','tenquille_lake','upper_squamish','wahleach_lake',
           'brookmere','callaghan','dickson_lake','disappointment_lake','dog_mountain','duffey_lake','gnawed_mountain',
           'grouse_mountain','hamilton_hill','highland_valley','klesilkwa','lightning_lake','mcgillivray_pass','nahatlatch',
           'orchid_lake','palisade_lake','shovelnose_mountain','stave_lake','sumallo_river_west','wahleach')
 
##'blackwall_peak_pillow','wahleach_lake',
sites <- c('blackwall_peak_pillow','chilliwack_river','spuzzum_creek','tenquille_lake','upper_squamish','wahleach_lake',
           'brookmere','callaghan','dickson_lake','disappointment_lake','dog_mountain','duffey_lake','gnawed_mountain',
           'great_bear','grouse_mountain','hamilton_hill','highland_valley','klesilkwa','lightning_lake','mcgillivray_pass',
           'nahatlatch','orchid_lake','palisade_lake','shovelnose_mountain','stave_lake','sumallo_river_west','wahleach')


##           'add_curd_mountain',
##           'add_scuzzy_mountain',
##           'add_russell_mountain',
##           'add_hemlock_valley',
##           'add_snowcap_peak',
##           'add_overseer_mountain',
##           'add_richmond',
##           'add_surrey',
##           'wahleach')

##Elevation data
prism.file <- '/storage/data/projects/rci/data/prism/van_whistler_prism_dem_elevations.nc'
prism.nc <- nc_open(prism.file)
dem <- ncvar_get(prism.nc,'elevation')
dem.lon <- ncvar_get(prism.nc,'lon')
dem.lat <- ncvar_get(prism.nc,'lat')
nc_close(prism.nc)
dem.diff <- (max(dem) - dem)/max(dem)
scale.elev <- 52 - dem.diff*12
slope.elev <- 0.88 - dem.diff*0.83
freq.elev <- 3.4 - dem.diff*0.6


##-----------------------------------------------------------------------
##Initial parameter values

param.names <- c('Lon','Lat','Elev','Dist','A','B','C','D')
site.params <- matrix(0,nrow=length(sites),ncol=length(param.names))

for (i in seq_along(sites)) {
   site <- sites[i]
   loc <- get_coordinates(site)
   dist <- (loc[1]+124)^2+(loc[2]-48.5)^2
   params <- get_cal_parameters(site,model)
   site.params[i,] <- c(loc,dist,params)
}
site.params[,5] <- site.params[,5]*-1

min.slopes <- site.params[,6] == 0.25
site.params[min.slopes,6] <- 0.35

##Grid for Kriging/TPS
##prism.file <- '/storage/data/projects/rci/data/winter_sports/PRISM/pr_monClim_PRISM_MODIS_GRID_198101-201012.nc' 
prism.file <- '/storage/data/projects/rci/data/winter_sports/PRISM/pr_monClim_PRISM_VAN_WHISTLER_198101-201012.nc' 
prism.nc <- nc_open(prism.file)
lon <- ncvar_get(prism.nc,'lon')
lat <- ncvar_get(prism.nc,'lat')
nc_close(prism.nc)

##if (model=='ERA5') {site.params[26,7] <- 2.2}

##site.params[site.params[,6] < 0.2,6]  <- 0.2

##---------------------------------------------------------
##Code for TPS
ras.xy <- cbind(c(matrix(dem.lon,nrow=nrow(dem),ncol=ncol(dem))),
                c(matrix(dem.lat,nrow=nrow(dem),ncol=ncol(dem),byrow=T)))

scale.tps <- Tps(x=data.frame(lon=site.params[,1],lat=site.params[,2],elev=site.params[,3]),
                 Y=site.params[,5])
scale.grid <- predict(scale.tps,x=cbind(ras.xy,c(dem)))
scale.matrix <- matrix(scale.grid,nrow=nrow(dem),ncol=ncol(dem))
##tps.ras <- list(x=dem.lon,y=dem.lat,z=tps.z)

slope.tps <- Tps(x=data.frame(lon=site.params[,1],lat=site.params[,2],elev=site.params[,3]),
                 Y=sqrt(site.params[,6]))
slope.grid <- predict(slope.tps,x=cbind(ras.xy,c(dem)))
slope.matrix <- (matrix(slope.grid,nrow=nrow(dem),ncol=ncol(dem)))^2 

freq.tps <- Tps(x=data.frame(lon=site.params[,1],lat=site.params[,2],elev=site.params[,3]),
                 Y=site.params[,7])
freq.grid <- predict(freq.tps,x=cbind(ras.xy,c(dem)))
freq.matrix <- matrix(freq.grid,nrow=nrow(dem),ncol=ncol(dem))


 write.dir <- '/storage/data/projects/rci/data/winter_sports/'
scale.file <- create_kriged_netcdf(var.name='scale',lon,lat,param=scale.matrix,
                                   model=model,write.dir=write.dir) 

slope.file <- create_kriged_netcdf(var.name='slope',lon,lat,param=slope.matrix,
                                   model=model,write.dir=write.dir) 
freq.file <- create_kriged_netcdf(var.name='freq',lon,lat,param=freq.matrix,
                                   model=model,write.dir=write.dir) 

browser()
##---------------------------------------------------------
##Code for ordinary kriging
if (1==0) {
  site.frame <- as.data.frame(site.params)
  names(site.frame) <- param.names
  coordinates(site.frame) <- ~ Lon + Lat

  scale.vgm <- variogram(log(A)~1,site.frame,cutoff=3)
  scale.fit <- fit.variogram(scale.vgm,model=vgm('Exp'))

  slope.vgm <- variogram(log(B)~1,site.frame,cutoff=3)
  slope.fit <- fit.variogram(slope.vgm,model=vgm('Exp'))

  freq.vgm <- variogram(log(C)~1,site.frame,cutoff=3)
  freq.fit <- fit.variogram(freq.vgm,model=vgm('Exp'))

  i <- expand.grid(seq(length(lon)), seq(length(lat)))
  names(i) <- c("xi", "yi")
  i.lon <- as.vector(lon[i$xi])
  i.lat <- as.vector(lat[i$yi])
  prism.coords <- data.frame(lon=i.lon, lat=i.lat, i)
  coordinates(prism.coords) <- c('lon','lat')

  scale.kriged <- krige(log(A)~1, site.frame, prism.coords, model=scale.fit)
  slope.kriged <- krige(log(B)~1, site.frame, prism.coords, model=slope.fit)
  freq.kriged <- krige(log(C)~1, site.frame, prism.coords, model=freq.fit)

  ##scale.matrix <- t(exp(matrix(scale.kriged@data$var1.pred,ncol=334,nrow=225,byrow=T)))
  ##slope.matrix <- t(exp(matrix(slope.kriged@data$var1.pred,ncol=334,nrow=225,byrow=T)))
  ##freq.matrix <- t(exp(matrix(freq.kriged@data$var1.pred,ncol=334,nrow=225,byrow=T)))

  scale.matrix <- t(exp(matrix(scale.kriged@data$var1.pred,ncol=480,nrow=323,byrow=T)))
  slope.matrix <- t(exp(matrix(slope.kriged@data$var1.pred,ncol=480,nrow=323,byrow=T)))
  freq.matrix <- t(exp(matrix(freq.kriged@data$var1.pred,ncol=480,nrow=323,byrow=T)))

  ##Plot parameters by Longitude, Latitude and Elevation

  scale.plus <- 0.3*scale.elev+0.7*scale.matrix
  slope.plus <- 0.3*slope.elev+0.7*slope.matrix
  freq.plus <- 0.3*freq.elev+0.7*freq.matrix

  write.dir <- '/storage/data/projects/rci/data/winter_sports/'
  ##scale.file <- create_kriged_netcdf(var.name='scale',lon,lat,param=scale.matrix,write.dir=write.dir) 
  ##slope.file <- create_kriged_netcdf(var.name='slope',lon,lat,param=slope.matrix,write.dir=write.dir) 
  ##freq.file <- create_kriged_netcdf(var.name='freq',lon,lat,param=freq.matrix,write.dir=write.dir) 

  par(mfrow=c(2,3),mar=c(3,3,1,1))
  image.plot(x=lon,y=lat,z=scale.matrix,
             xlim=c(-123.7,-120.65),ylim=c(48.85,50.8))
  points(x=site.params[,1],y=site.params[,2],pch=18)
  contour(x=lon,y=lat,z=scale.matrix,nlevels=10,add=T)

  image.plot(x=lon,y=lat,z=slope.matrix,
           xlim=c(-123.7,-120.65),ylim=c(48.85,50.8))
  points(x=site.params[,1],y=site.params[,2],pch=18)
  contour(x=lon,y=lat,z=slope.matrix,nlevels=10,add=T)

  image.plot(x=lon,y=lat,z=freq.matrix,
             xlim=c(-123.7,-120.65),ylim=c(48.85,50.8))
  points(x=site.params[,1],y=site.params[,2],pch=18)
  contour(x=lon,y=lat,z=freq.matrix,nlevels=10,add=T)

  image.plot(x=lon,y=lat,z=scale.plus,
             xlim=c(-123.7,-120.65),ylim=c(48.85,50.8))
  points(x=site.params[,1],y=site.params[,2],pch=18)
  contour(x=lon,y=lat,z=scale.plus,nlevels=10,add=T)

  image.plot(x=lon,y=lat,z=slope.plus,
           xlim=c(-123.7,-120.65),ylim=c(48.85,50.8))
  points(x=site.params[,1],y=site.params[,2],pch=18)
  contour(x=lon,y=lat,z=slope.plus,nlevels=10,add=T)

  image.plot(x=lon,y=lat,z=freq.plus,
           xlim=c(-123.7,-120.65),ylim=c(48.85,50.8))
  points(x=site.params[,1],y=site.params[,2],pch=18)
  contour(x=lon,y=lat,z=freq.plus,nlevels=10,add=T)
}



browser()

##Longitude

lons <- as.numeric(site.params[,1])
lats <- as.numeric(site.params[,2])
elevs <- as.numeric(site.params[,3])
if (1==1) {

x11()
par(mfrow=c(3,1))
plot(lons[order(lons)],as.numeric(site.params[order(lons),4]),type='l',lwd=2,main='SCF')
points(lons[order(lons)],as.numeric(site.params[order(lons),4]),pch=18,cex=2,main='SCF')
plot(lons[order(lons)],as.numeric(site.params[order(lons),5]),type='l',lwd=2,main='MFMax')
points(lons[order(lons)],as.numeric(site.params[order(lons),5]),pch=18,cex=2)
plot(lons[order(lons)],as.numeric(site.params[order(lons),6]),type='l',lwd=2,main='TIPM')
points(lons[order(lons)],as.numeric(site.params[order(lons),6]),pch=18,cex=2)

x11()
par(mfrow=c(3,1))
plot(lats[order(lats)],as.numeric(site.params[order(lats),4]),type='l',lwd=2,main='SCF')
points(lats[order(lats)],as.numeric(site.params[order(lats),4]),pch=18,cex=2)
plot(lats[order(lats)],as.numeric(site.params[order(lats),5]),type='l',lwd=2,main='MFMax')
points(lats[order(lats)],as.numeric(site.params[order(lats),5]),pch=18,cex=2)
plot(lats[order(lats)],as.numeric(site.params[order(lats),6]),type='l',lwd=2,main='TIPM')
points(lats[order(lats)],as.numeric(site.params[order(lats),6]),pch=18,cex=2)

x11()
par(mfrow=c(3,1))
plot(elevs[order(elevs)],as.numeric(site.params[order(elevs),4]),type='l',lwd=2,main='SCF')
points(elevs[order(elevs)],as.numeric(site.params[order(elevs),4]),pch=18,cex=2)
plot(elevs[order(elevs)],as.numeric(site.params[order(elevs),5]),type='l',lwd=2,main='MFMax')
points(elevs[order(elevs)],as.numeric(site.params[order(elevs),5]),pch=18,cex=2)
plot(elevs[order(elevs)],as.numeric(site.params[order(elevs),6]),type='l',lwd=2,main='TIPM')
points(elevs[order(elevs)],as.numeric(site.params[order(elevs),6]),pch=18,cex=2)


}