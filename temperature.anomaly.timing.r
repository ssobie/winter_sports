##Find the years (seasons?) when 1,2,3 degrees of warming have occurred 
##in the Van-Whistler region

library(raster)
library(zoo)

##---------------------------------------------------------
##Create time series
create_time_series <- function(tasmax.file,tasmin.file,map.extent,tmp.dir) {
   tasmax <- brick(paste0(tmp.dir,tasmax.file))
   tasmin <- brick(paste0(tmp.dir,tasmin.file))
   tas <- (tasmax + tasmin) / 2
   tas.crop <- crop(tas,map.extent)
   tas.series <- as.numeric(cellStats(tas.crop,mean,na.rm=T))
   dates <- tasmax@z$Date
   tas.zoo <- zoo(tas.series,dates)
   tas.roll <- rollmean(tas.zoo,31)
   return(list(roll=tas.roll,raw=tas.zoo))
}

##---------------------------------------------------------

find_anomaly_timing <- function(tas.series,base,anoms) {

   tas.raw <- tas.series$raw
   raw.dates <- index(tas.raw)
   raw.years <- format(raw.dates,'%Y')
   raw.series <- as.numeric(tas.raw)

   rst <- head(grep(base[1],raw.years),1)
   ren <- tail(grep(base[2],raw.years),2)
   raw.baseline <- mean(raw.series[rst:ren])
   raw.anomalies <- raw.series - raw.baseline

   tas.roll <- tas.series$roll
   dates <- index(tas.roll)
   years <- format(dates,'%Y')
   series <- as.numeric(tas.roll)

   bst <- head(grep(base[1],years),1)
   ben <- tail(grep(base[2],years),2)
   baseline <- mean(series[bst:ben])
   anomalies <- series - baseline

   anom.years <- rep(0,length(anoms))
   anom.values <- rep(0,length(anoms))

   for (i in seq_along(anoms)) {
      ix <- which.min(abs(anoms[i] - anomalies))
      anom.years[i] <- years[ix]

      rx <- which(years[ix] == raw.years)      
      anom.values[i] <- mean(raw.anomalies[(rx-15):(rx+15)])
      if (is.na(anom.values[i])) { browser()}
   }   
   return(list(years=anom.years,values=anom.values,anoms=raw.anomalies,rolled=anomalies))
}


##---------------------------------------------------------

data.dir <- "/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/"

tmp.dir <- "/local_temp/ssobie/snow_temp_timing/"

if (!file.exists(tmp.dir)) {
   dir.create(tmp.dir,recursive=T)
}

gcm.list <- c('ACCESS1-0','CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G',
              'HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MPI-ESM-LR','MRI-CGCM3')
map.extent <- extent(c(-123.7,-120.65,48.9,50.8))
anoms <- c(1,2,3,4)
base <- c(1981,2010)

raw.anomalies <- vector(mode='list',length=length(gcm.list))
rolled.anomalies <- vector(mode='list',length=length(gcm.list))

anomaly.years <- matrix(0,nrow=length(gcm.list),ncol=length(anoms))
anomaly.values <- matrix(0,nrow=length(gcm.list),ncol=length(anoms))

for (j in seq_along(gcm.list)) {

   gcm <- gcm.list[j]
   print(gcm)
   gcm.dir <- paste0(data.dir,gcm,'/tas_climatology/')
   tasmax.file <- list.files(path=gcm.dir,pattern=paste0("tasmax_annual_mean_gcm_prism_",gcm))
   tasmin.file <- list.files(path=gcm.dir,pattern=paste0("tasmin_annual_mean_gcm_prism_",gcm))

   file.copy(from=paste0(gcm.dir,tasmax.file),to=tmp.dir,overwrite=TRUE)
   file.copy(from=paste0(gcm.dir,tasmin.file),to=tmp.dir,overwrite=TRUE)

   tas.roll <- create_time_series(tasmax.file,tasmin.file,map.extent,tmp.dir)
   tas.timing <- find_anomaly_timing(tas.roll,base,anoms)
   anomaly.years[j,] <- tas.timing$years
   anomaly.values[j,] <- tas.timing$values
   raw.anomalies[[j]] <- tas.timing$anoms
   rolled.anomalies[[j]] <- tas.timing$rolled
   
   file.remove(paste0(tmp.dir,tasmax.file))
   file.remove(paste0(tmp.dir,tasmin.file))

}
rv <- rbind(c('Model','1Degrees','2Degrees','3Degrees','4Degrees'),cbind(gcm.list,anomaly.years))
write.table(rv,file='/storage/data/projects/rci/data/winter_sports/temperature.anomaly.years.csv',
            sep=',',quote=FALSE,row.name=FALSE,col.name=F)
 
##plot(rolled.anomalies[[1]],type='l',ylim=c(-2,8))
##for (j in seq_along(gcm.list)) {
##   lines(rolled.anomalies[[j]])
##}