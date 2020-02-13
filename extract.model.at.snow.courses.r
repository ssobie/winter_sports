##Script to pull out the time series of pr, tasmax, tamsin from the driving models 
library(ncdf4)
library(PCICt)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

get.coordinates <- function(site) {

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
                      spuzzum_creek=c(-121.686,49.74,1197),
                      chilliwack_river=c(-121.71667,49.0333,1600),
                      upper_squamish=c(-123.4333,50.1500,1340),
                      wahleach_lake=c(-121.5833,49.2333,1400),
                      tenquille_lake=c(-122.9333,50.5333,1680))

  rv <- coordinates[[site]]
  return(rv)
}

get.800m.data <- function(site, pr.nc,tasmax.nc,tasmin.nc) {

  coords <- get.coordinates(site)
  plot.title <- site
  ##
  lon <- ncvar_get(pr.nc,'lon')
  lat <- ncvar_get(pr.nc,'lat')

  lon.bnds <- coords[1]
  lat.bnds <- coords[2]
  elev <- coords[3]

  lon.ix <- which.min(abs(lon-lon.bnds)) ##196 ####192, 61
  lat.ix <- which.min(abs(lat-lat.bnds)) ##55 ##
  print(lon.ix)
  print(lon[lon.ix])
  print(lat.ix)
  print(lat[lat.ix])


  pr.data <- ncvar_get(pr.nc,'pr',start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
  tasmax.data <- ncvar_get(tasmax.nc,'tasmax',start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
  tasmin.data <- ncvar_get(tasmin.nc,'tasmin',start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
  flag <- which(tasmax.data<tasmin.data)
  tasmax.data[flag] <- tasmin.data[flag]+1.1
  tas.data <- (tasmax.data+tasmin.data)/2

  rv <- list(pr=pr.data,
             tasmax=tasmax.data,
             tasmin=tasmin.data,
             tas=tas.data)

  return(rv)
}

  sites <- c('callaghan',
             'orchid_lake',
             'palisade_lake',
             'grouse_mountain',
             'dog_mountain',
             'stave_lake',
             'nahatlatch',
             'wahleach',
             'klesilkwa',
             'lightning_lake',
             'brookmere',
             'shovelnose_mountain',
             'hamilton_hill',
             'spuzzum_creek',
             'chilliwack_river',
             'upper_squamish',
             'wahleach_lake',             
             'tenquille_lake')
sites <- 'grouse_station'

base.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/'
###base.dir <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/'

model <- 'PNWNAmet'

###pr.file <- paste0('pr_gcm_prism_',model,'_19790101-20181031.nc')
###tasmax.file <- paste0('tasmax_gcm_prism_',model,'_19790101-20181031.nc')
###tasmin.file <- paste0('tasmin_gcm_prism_',model,'_19790101-20181031.nc')
pr.file <- paste0('pr_gcm_prism_allBC_TPS_1945-2012.nc')
tasmax.file <- paste0('tasmax_gcm_prism_allBC_TPS_1945-2012.nc')
tasmin.file <- paste0('tasmin_gcm_prism_allBC_TPS_1945-2012.nc')

pr.nc <- nc_open(paste(base.dir,model,'/',pr.file,sep=''))
tasmax.nc <- nc_open(paste(base.dir,model,'/',tasmax.file,sep=''))
tasmin.nc <- nc_open(paste(base.dir,model,'/',tasmin.file,sep=''))

dates <- netcdf.calendar(pr.nc)

  for (site in sites) {
    print(site)
    clim.data <- get.800m.data(site, pr.nc,tasmax.nc,tasmin.nc)
    output <- cbind(as.character(dates),round(cbind(clim.data$pr,clim.data$tasmax,
                               clim.data$tasmin,clim.data$tas),2))
    output <- rbind(c('Dates','Pr','Tasmax','Tasmin','Tas'),output)
    write.table(output,file=paste('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/',site,'_',model,'_800m_data.csv',sep=''),
                sep=',',row.name=FALSE,col.name=FALSE,quote=FALSE)

  }
nc_close(pr.nc)
nc_close(tasmax.nc)
nc_close(tasmin.nc)