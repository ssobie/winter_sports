##Script to compare snowdepth, SWE and snow cover between modelled values, 
##snow course sites and MODIS snow cover
library(ncdf4)
library(PCICt)

source('/storage/home/ssobie/code/repos/winter_sports/test.snow.model.r',chdir=T)
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
                      spuzzum_creek=c(-121.686,49.674,1197))
  rv <- coordinates[[site]]
  return(rv)
}


get.tas.data <- function(model,suffix,coords) {
             
    model.tx.nc <- nc_open(paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/',model,'/tasmax_',suffix))
    lon <- ncvar_get(model.tx.nc,'lon')
    lat <- ncvar_get(model.tx.nc,'lat')
    lon.ix <- which.min(abs(coords[1]-lon))
    lat.ix <- which.min(abs(coords[2]-lat))
    print('Location indices')
    print(c(lon.ix,lat.ix))
        
    model.tx <- ncvar_get(model.tx.nc,'tasmax',start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
    model.time <- netcdf.calendar(model.tx.nc)
    nc_close(model.tx.nc)
    model.tn.nc <- nc_open(paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/',model,'/tasmin_',suffix))
    model.tn <- ncvar_get(model.tn.nc,'tasmin',start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
    nc_close(model.tn.nc)

    model.tas <- (model.tx+model.tn)/2

    rv <- list(tas=model.tas,
               time=as.Date(as.character(model.time)))
    return(rv)
}

site <- 'spuzzum_creek'
model <- 'NCEP2'

coords <- get.coordinates(site)
lat.bnds <- coords[2]
elev <- coords[3]


##Loop over sites

    ncep2.raw <- get.tas.data('NCEP2','day_QDM_NCEP2_1979-2016.nc',coords)
    ncep2.anoms <- get.tas.data('NCEP2','anoms_QDM_NCEP2_1979-2016.nc',coords)
    ncep2.anoms.interp <- get.tas.data('NCEP2','anoms_interp_NCEP2_1979-2016.nc',coords)
    ncep2.gcm.prism <- get.tas.data('NCEP2','gcm_prism_NCEP2_1979-2016.nc',coords)
    

    ##Snow Pillow Data
    pillow.file <- paste('/storage/data/projects/rci/data/assessments/snow_model/snow_pillow/',site,'_asp.csv',sep='')
    pillow.data <- read.csv(pillow.file,header=T,as.is=T)
    pillow.dates <- format(as.Date(pillow.data[,2]),'%Y-%m-%d')
    pillow.tasmax <- pillow.data[,3]
    pillow.tasmin <- pillow.data[,5]
    pillow.tas <- (pillow.tasmax + pillow.tasmin)/2
    pillow.precip <- pillow.data[,7]##mm
    sb <- 1:1000

    date.subset <- format(ncep2.raw$time,'%Y-%m-%d') %in% pillow.dates[sb]
    
    par(mfrow=c(2,2))    

    plot(ncep2.raw$time[date.subset],round(ncep2.raw$tas[date.subset],1),type='l',lwd=3,col='red',main='NCEP2 Raw TAS (C)',cex.axis=1.5)
    lines(as.Date(pillow.dates)[sb],pillow.tas[sb],lwd=3,col='orange')
    abline(h=0)

    plot(ncep2.anoms$time[date.subset],round(ncep2.anoms$tas[date.subset],1),type='l',lwd=3,col='red',main='NCEP ANOMS TAS (C)',cex.axis=1.5)
    ##lines(as.Date(pillow.dates)[sb],pillow.tas[sb],lwd=3,col='orange')
    abline(h=0)

    plot(ncep2.anoms.interp$time[date.subset],round(ncep2.anoms.interp$tas[date.subset],1),type='l',lwd=3,col='red',
         main='NCEP ANOMS INTERP TAS (C)',cex.axis=1.5)
    ##lines(as.Date(pillow.dates)[sb],pillow.tas[sb],lwd=3,col='orange')
    abline(h=0)

    plot(ncep2.gcm.prism$time[date.subset],round(ncep2.gcm.prism$tas[date.subset],1),type='l',lwd=3,col='red',
         main='NCEP GCM-PRISM TAS (C)',cex.axis=1.5)
    lines(as.Date(pillow.dates)[sb],pillow.tas[sb],lwd=3,col='orange')
    abline(h=0)

