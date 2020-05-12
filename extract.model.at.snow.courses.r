##Script to pull out the time series of pr, tasmax, tamsin from the driving models 
library(ncdf4)
library(PCICt)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)
source('/storage/home/ssobie/code/repos/winter_sports/site.coordinates.r',chdir=TRUE)

##----------------------------------------------------------------------------

get_800m_data <- function(site, pr.nc,tasmax.nc,tasmin.nc) {

  coords <- get_coordinates(site)
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

##----------------------------------------------------------------------------
get_snow_model_data <- function(site, swe.nc, snowdepth.nc) {

  coords <- get_coordinates(site)
  plot.title <- site
  ##
  lon <- ncvar_get(swe.nc,'lon')
  lat <- ncvar_get(swe.nc,'lat')

  lon.bnds <- coords[1]
  lat.bnds <- coords[2]
  elev <- coords[3]

  lon.ix <- which.min(abs(lon-lon.bnds)) ##196 ####192, 61
  lat.ix <- which.min(abs(lat-lat.bnds)) ##55 ##
  print(lon.ix)
  print(lon[lon.ix])
  print(lat.ix)
  print(lat[lat.ix])


  swe.data <- ncvar_get(swe.nc,'swe',start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
  depth.data <- ncvar_get(snowdepth.nc,'snowdepth',start=c(lon.ix,lat.ix,1),count=c(1,1,-1))

  rv <- list(swe=swe.data,
             depth=depth.data)
  return(rv)
}


##----------------------------------------------------------------------------

  sites <- c('blackwall_peak_course',
             'blackwall_peak_pillow',
             'boston_bar_lower',
             'boston_bar_upper',
             'brookmere', 
             'burwell_lake',
             'callaghan',
             'chapman_creek',
             'chilliwack_river', 
             'cornwall_hills',
             'diamond_head',
             'dickson_lake',
             'disappointment_lake',
             'dog_mountain',
             'duffey_lake',
             'edwards_lake',
             'garibaldi_lake',
             'gnawed_mountain',
             'great_bear',
             'grouse_mountain',
             'hamilton_hill',
             'highland_valley',
             'hollyburn',
             'hope',
             'klesilkwa',
             'lightning_lake',             
             'loch_lomond',
             'lytton',
             'mcgillivray_pass',
             'mount_seymour',
             'nahatlatch',
             'new_tashme',                                     
             'orchid_lake',
             'ottomite',
             'palisade_lake',
             'pavilion_mountain',
             'shalalth',
             'shovelnose_mountain',
             'spuzzum_creek',             
             'stave_lake',
             'sumallo_river',
             'sumallo_river_west',
             'tenquille_lake',
             'tenquille_course',
             'upper_squamish',
             'wahleach',
             'wahleach_lake',             
             'whistler_mountain',         
             'wolverine_creek')



##Extract the precip and temperature series
if (1==0) {

base.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/'
###base.dir <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/'

model <- 'PNWNAmet'

###pr.file <- paste0('pr_gcm_prism_',model,'_BC_19800101-20181231.nc')
###tasmax.file <- paste0('tasmax_gcm_prism_',model,'_BC_19800101-20181231.nc')
###tasmin.file <- paste0('tasmin_gcm_prism_',model,'_BC_19800101-20181231.nc')

pr.file <- paste0('pr_gcm_prism_allBC_TPS_1945-2012.nc')
tasmax.file <- paste0('tasmax_gcm_prism_allBC_TPS_1945-2012.nc')
tasmin.file <- paste0('tasmin_gcm_prism_allBC_TPS_1945-2012.nc')

pr.nc <- nc_open(paste(base.dir,model,'/',pr.file,sep=''))
tasmax.nc <- nc_open(paste(base.dir,model,'/',tasmax.file,sep=''))
tasmin.nc <- nc_open(paste(base.dir,model,'/',tasmin.file,sep=''))

dates <- netcdf.calendar(pr.nc)

  for (site in sites) {
    print(site)
    clim.data <- get_800m_data(site, pr.nc,tasmax.nc,tasmin.nc)
    output <- cbind(as.character(dates),round(cbind(clim.data$pr,clim.data$tasmax,
                               clim.data$tasmin,clim.data$tas),2))
    output <- rbind(c('Dates','Pr','Tasmax','Tasmin','Tas'),output)
    write.table(output,file=paste('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/',site,'_',model,'_800m_data.csv',sep=''),
                sep=',',row.name=FALSE,col.name=FALSE,quote=FALSE)

  }
nc_close(pr.nc)
nc_close(tasmax.nc)
nc_close(tasmin.nc)

browser()

}



##-----------------------------------------------------------------------------

##base.dir <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/calibrated_PNWNAmet_prism_with_elevation/'
##model <- 'PNWNAmet'
##type <- 'PRISM_with_elevation'
##swe.file <- paste0('swe_BCCAQ2-PRISM_',model,'_1945-2012.nc')
##snowdepth.file <- paste0('snowdepth_BCCAQ2-PRISM_',model,'_1945-2012.nc')

##------------------------------------

##base.dir <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/calibrated_PNWNAmet_prism_tps/'
base.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/snow_model/'
model <- 'MRI-CGCM3'
reanalysis <- 'PNWNAmet'
type <- 'PRISM_TPS'

read.dir <- paste0(base.dir,'calibrated_',model,'_',reanalysis,'_prism_tps/')
swe.file <- paste0('swe_BCCAQ2-PRISM_',model,'_',reanalysis,'_19500101-21001231.nc')
snowdepth.file <- paste0('snowdepth_BCCAQ2-PRISM_',model,'_',reanalysis,'_19500101-21001231.nc')

swe.nc <- nc_open(paste(read.dir,swe.file,sep=''))
snowdepth.nc <- nc_open(paste(read.dir,snowdepth.file,sep=''))

dates <- netcdf.calendar(swe.nc)

  for (site in sites) {
    print(site)
    snow.data <- get_snow_model_data(site, swe.nc,snowdepth.nc)
    output <- cbind(as.character(dates),round(cbind(snow.data$swe,snow.data$depth),2))
                               
    output <- rbind(c('Dates','SWE','Snowdepth'),output)

    write.table(output,
                file=paste('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow_series/',
                           site,'_',model,'_',reanalysis,'_',type,'_snow_model_data.csv',sep=''),
                           sep=',',row.name=FALSE,col.name=FALSE,quote=FALSE)
  }
nc_close(swe.nc)
nc_close(snowdepth.nc)

