##Snow-17 accumulation and ablation model.

library(ncdf4)
library(sp)
library(gstat)
library(fields)

##-------------------------------------------------------------------------------
##Create netcdfs for kriged parameters

create_param_netcdf <- function(var.name,lon,lat,param,write.dir) {

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

  write.clim.name <- paste0(write.dir,var.name,'_hyper_snow_calibrated_parameter_with_elevation.nc')

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
##Grid
prism.file <- '/storage/data/projects/rci/data/prism/bc_prism_MODIS_grid_dem_elevations.nc'
prism.nc <- nc_open(prism.file)
lon <- ncvar_get(prism.nc,'lon')
lat <- ncvar_get(prism.nc,'lat')
dem <- ncvar_get(prism.nc,'elevation')
dem.raster <- list(x=lon,y=lat,z=dem)
dem.diff <- (max(dem) - dem)/max(dem)

##elev.raster <- list(x=lon,y=lat,z=scale.elev-46)

scale.file <- '/storage/data/projects/rci/data/winter_sports/scale_hyper_snow_calibrated_parameter.nc'
scale.nc <- nc_open(scale.file)
scale.data <- ncvar_get(scale.nc,'scale')
scale.elev <- 52 - dem.diff*12
scale.raster <- list(x=lon,y=lat,z=scale.data)
nc_close(scale.nc)

##scale.plus <- list(x=lon,y=lat,z=0.4*scale.elev+0.6*scale.data)
write.dir <- '/storage/data/projects/rci/data/winter_sports/'
scale.plus <- 0.3*scale.elev+0.7*scale.data
scale.file <- create_param_netcdf(var.name='scale',lon,lat,param=scale.plus,write.dir=write.dir) 

slope.file <- '/storage/data/projects/rci/data/winter_sports/slope_hyper_snow_calibrated_parameter.nc'
slope.nc <- nc_open(slope.file)
slope.data <- ncvar_get(slope.nc,'slope')
slope.elev <- 0.88 - dem.diff*0.83
slope.plus <- 0.3*slope.elev +0.7*slope.data
slope.file <- create_param_netcdf(var.name='slope',lon,lat,param=slope.plus,write.dir=write.dir) 
nc_close(slope.nc)

freq.file <- '/storage/data/projects/rci/data/winter_sports/freq_hyper_snow_calibrated_parameter.nc'
freq.nc <- nc_open(freq.file)
freq.data <- ncvar_get(freq.nc,'freq')
freq.elev <- 3.4 - dem.diff*0.6
freq.plus <- 0.3*freq.elev +0.7*freq.data
freq.file <- create_param_netcdf(var.name='freq',lon,lat,param=freq.plus,write.dir=write.dir) 
nc_close(freq.nc)



image.plot(scale.raster)
x11()

image.plot(scale.plus)

radius <- dem*0
for (i in 1:length(lon)) {
   for (j in 1:length(lat)) {
      radius[i,j] <- sqrt( ((lon[i]-min(lon)+4)^2)/6 + ((lat[j]-min(lat))^2))
   }
}

scaled.radius <- radius/max(radius)


cosine <- cos(scaled.radius*pi/2)
scale.cosine <- 40 + cosine*12
#image.plot(scale.cosine)

mean.scale <- (scale.elev+scale.cosine)/2

space <- list(x=lon,y=lat,z=mean.scale)
#image.plot(space)

#i <- expand.grid(seq(length(lon)), seq(length(lat)))
#names(i) <- c("xi", "yi")
#i.lon <- as.vector(lon[i$xi])
#i.lat <- as.vector(lat[i$yi])
#prism.coords <- data.frame(lon=i.lon, lat=i.lat, i)
#coordinates(prism.coords) <- c('lon','lat')

nc_close(prism.nc)



##scale.matrix <- t(exp(matrix(scale.kriged@data$var1.pred,ncol=334,nrow=225,byrow=T)))
##slope.matrix <- t(exp(matrix(slope.kriged@data$var1.pred,ncol=334,nrow=225,byrow=T)))
##freq.matrix <- t(exp(matrix(freq.kriged@data$var1.pred,ncol=334,nrow=225,byrow=T)))

##write.dir <- '/storage/data/projects/rci/data/winter_sports/'
##scale.file <- create_kriged_netcdf(var.name='scale',lon,lat,param=scale.matrix,write.dir=write.dir) 
##slope.file <- create_kriged_netcdf(var.name='slope',lon,lat,param=slope.matrix,write.dir=write.dir) 
##freq.file <- create_kriged_netcdf(var.name='freq',lon,lat,param=freq.matrix,write.dir=write.dir) 




##Plot parameters by Longitude, Latitude and Elevation

#image.plot(x=lon,y=lat,z=scale.matrix,
#           xlim=c(-123.55,-120.6),ylim=c(48.9,50.7))
#contour(x=lon,y=lat,z=scale.matrix,nlevels=10,add=T)

