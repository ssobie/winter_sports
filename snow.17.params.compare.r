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

##---------------------------------------------------------------------
save.dir <- '/storage/data/projects/rci/data/winter_sports/obs/snow17/'
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
           ##SCF, MFMax,  TIPM
sub.par <- c(0.75, 0.75,  0.05)

param.names <- c('Lon','Lat','Elev','Pr','Snow','TIPM')
site.params <- matrix(0,nrow=length(sites),ncol=length(param.names))

for (i in seq_along(sites)) {
   site <- sites[i]
   loc <- get_coordinates(site)
   load(paste0(save.dir,site,'_snow17_optim_3_parameter_fit_bounded_pnwnamet.RData'))
   site.params[i,] <- c(loc,optim.result$par)
   rm(optim.result)
}

site.frame <- as.data.frame(site.params)
names(site.frame) <- param.names
coordinates(site.frame) <- ~ Lon + Lat

pr.vgm <- variogram(log(Pr)~1,site.frame)
pr.fit <- fit.variogram(pr.vgm,model=vgm('Bes'))

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

pr.kriged <- krige(log(Pr)~1, site.frame, prism.coords, model=pr.fit)

##Plot parameters by Longitude, Latitude and Elevation

image.plot(x=lon,y=lat,z=t(exp(matrix(pr.kriged@data$var1.pred,ncol=334,nrow=225,byrow=T))),
           xlim=c(-123.55,-120.6),ylim=c(48.9,50.7))
contour(x=lon,y=lat,z=t(exp(matrix(pr.kriged@data$var1.pred,ncol=334,nrow=225,byrow=T))),nlevels=10,add=T)
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