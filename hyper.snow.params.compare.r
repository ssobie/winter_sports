##Snow-17 accumulation and ablation model.

library(ncdf4)
library(sp)
library(gstat)
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

create_kriged_netcdf <- function(var.name,lon,lat,param,write.dir) {

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

  write.clim.name <- paste0(write.dir,var.name,'_hyper_snow_calibrated_parameter.nc')

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



get_coordinates <- function(site) {

  coordinates <- list(callaghan=c(-123.1036,50.1383278,1009),
                      orchid_lake=c(-123.0519638,49.53678,1178),
                      palisade_lake=c(-123.0321944,49.454433,898),
                      grouse_mountain=c(-123.0774472,49.383655,1126),
                      dog_mountain=c(-122.96255,49.37251944,1007),
                      dickson_lake=c(-122.06984166,49.3168194,1147),
                      stave_lake=c(-122.315805,49.58030277,1211),
                      nahatlatch=c(-122.059261,49.825866,1530),
                      wahleach=c(-121.57945,49.2298694,1395),
                      klesilkwa=c(-121.3086527,49.129438,610),
                      lightning_lake=c(-120.850205,49.044788,1254),
                      brookmere=c(-120.87397,49.815027,994),
                      shovelnose_mountain=c(-120.864175,49.8546305,1456),
                      hamilton_hill=c(-120.7955805,49.4988027,1477),
                      spuzzum_creek=c(-121.686,49.674,1197),
                      chilliwack_river=c(-121.71667,49.0333,1600),
                      upper_squamish=c(-123.4333,50.1500,1340),
                      wahleach_lake=c(-121.5833,49.2333,1400),
                      tenquille_lake=c(-122.9333,50.5333,1680))

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
                    brookmere=c(40,0.05,2.7,1.0209),
                    shovelnose_mountain=c(40.15,0.075,2.7,1.0209),
                    lightning_lake=c(40.15,0.05,2.75,1.0209),
                    hamilton_hill=c(40.45,0.1,2.85,1.0209),
                    upper_squamish=c(51.7,0.575,3.7,1.0209),
                    tenquille_lake=c(43.45,0.5,2.25,1.0209),
                    chilliwack_river=c(51.8,0.8,4.0,1.0209),
                    spuzzum_creek=c(51.8,0.875,4.1,1.0209))
   rv <- site.params[[site]]
   return(rv)
}

##---------------------------------------------------------------------

sites <- c('shovelnose_mountain',
           'brookmere',
           'lightning_lake',
           'callaghan',
           'orchid_lake',
           'palisade_lake',
           'grouse_mountain',
           'dog_mountain',
           'stave_lake',
           'nahatlatch',
           'wahleach',
           'klesilkwa',
           'hamilton_hill',
           'upper_squamish',
           'spuzzum_creek',
           'chilliwack_river',
           'tenquille_lake')



##-----------------------------------------------------------------------
##Initial parameter values

param.names <- c('Lon','Lat','Elev','A','B','C','D')
site.params <- matrix(0,nrow=length(sites),ncol=length(param.names))

for (i in seq_along(sites)) {
   site <- sites[i]
   loc <- get_coordinates(site)
   params <- get_parameters(site)
   site.params[i,] <- c(loc,params)
}

site.frame <- as.data.frame(site.params)
names(site.frame) <- param.names
coordinates(site.frame) <- ~ Lon + Lat

scale.vgm <- variogram(log(A)~1,site.frame,cutoff=3)
scale.fit <- fit.variogram(scale.vgm,model=vgm('Exp'))

slope.vgm <- variogram(log(B)~1,site.frame,cutoff=2)
slope.fit <- fit.variogram(slope.vgm,model=vgm('Exp'))

freq.vgm <- variogram(log(C)~1,site.frame,cutoff=4)
freq.fit <- fit.variogram(freq.vgm,model=vgm('Exp'))


##Grid for Kriging
prism.file <- '/storage/data/projects/rci/data/winter_sports/PRISM/pr_monClim_PRISM_MODIS_GRID_198101-201012.nc' 
prism.nc <- nc_open(prism.file)
lon <- ncvar_get(prism.nc,'lon')
lat <- ncvar_get(prism.nc,'lat')

i <- expand.grid(seq(length(lon)), seq(length(lat)))
names(i) <- c("xi", "yi")
i.lon <- as.vector(lon[i$xi])
i.lat <- as.vector(lat[i$yi])
prism.coords <- data.frame(lon=i.lon, lat=i.lat, i)
coordinates(prism.coords) <- c('lon','lat')

nc_close(prism.nc)

scale.kriged <- krige(log(A)~1, site.frame, prism.coords, model=scale.fit)
slope.kriged <- krige(log(B)~1, site.frame, prism.coords, model=slope.fit)
freq.kriged <- krige(log(C)~1, site.frame, prism.coords, model=freq.fit)

scale.matrix <- t(exp(matrix(scale.kriged@data$var1.pred,ncol=334,nrow=225,byrow=T)))
slope.matrix <- t(exp(matrix(slope.kriged@data$var1.pred,ncol=334,nrow=225,byrow=T)))
freq.matrix <- t(exp(matrix(freq.kriged@data$var1.pred,ncol=334,nrow=225,byrow=T)))

write.dir <- '/storage/data/projects/rci/data/winter_sports/'

scale.file <- create_kriged_netcdf(var.name='scale',lon,lat,param=scale.matrix,write.dir=write.dir) 
slope.file <- create_kriged_netcdf(var.name='slope',lon,lat,param=slope.matrix,write.dir=write.dir) 
freq.file <- create_kriged_netcdf(var.name='freq',lon,lat,param=freq.matrix,write.dir=write.dir) 

##Plot parameters by Longitude, Latitude and Elevation

image.plot(x=lon,y=lat,z=scale.matrix,
           xlim=c(-123.55,-120.6),ylim=c(48.9,50.7))
contour(x=lon,y=lat,z=scale.matrix,nlevels=10,add=T)
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