##Script to compare snowdepth, SWE and snow cover between modelled values, 
##snow course sites and MODIS snow cover

##Rprof('profile.snow.model.out')

library(ncdf4)
library(plotrix)
library(zoo)

source('/storage/home/ssobie/code/repos/winter_sports/snow.model.functional.r',chdir=T)
##source('/storage/home/ssobie/code/repos/winter_sports/test.snow.model.r',chdir=T)

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
                      spuzzum_creek=c(-121.686,49.674,1197),
                      chilliwack_river=c(-121.71667,49.0333,1600),
                      upper_squamish=c(-123.4333,50.1500,1340),
                      wahleach_lake=c(-121.5833,49.2333,1400),
                      tenquille_lake=c(-122.9333,50.5333,1680))

  rv <- coordinates[[site]]
  return(rv)
}

hyper.snow <- function(pr,tasmax,tasmin,coeffs) {
        
        tas <- (tasmax+tasmin)/2                # degrees C
        frac <- coeffs$a*(tanh(coeffs$b*(tas-coeffs$c))-coeffs$d)
        sample <- runif(length(tas),min=0,max=100)
        test <- sample > frac        
        high.temp <- tas > 12
        test[high.temp] <- TRUE
        snow.type <- rep(TRUE,length(tas))
        snow.type[test] <- FALSE

        NewSnowWatEq <- pr*0.001
        NewSnowWatEq[!snow.type] <- 0
        R_m <- pr*0.001
        R_m[snow.type] <- 0
        rv <- list(snow=NewSnowWatEq,
                   rain=R_m)
        return(rv)
}


##Slope and Aspect Values
as.dir <- '/storage/data/projects/rci/data/prism/'
slopes.nc <- nc_open(paste0(as.dir,'prism_slopes.nc'))
bc.slopes <- ncvar_get(slopes.nc,'Band1')/90*pi/2
bc.lon <- ncvar_get(slopes.nc,'lon')
bc.lat <- ncvar_get(slopes.nc,'lat')
nc_close(slopes.nc)

aspects.nc <- nc_open(paste0(as.dir,'prism_aspects.nc')) 
bc.aspects <- ncvar_get(aspects.nc,'Band1')/360*2*pi
nc_close(aspects.nc)

slen <- 1001
save.dir <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sims/'

ncep2.file <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/grouse_mountain_NCEP2_800m_data.csv'
ncep2.data <- read.csv(ncep2.file,header=T,as.is=T)

ncep2.swe.sims <- matrix(0,nrow=dim(ncep2.data)[1],ncol=slen)
ncep2.snow.sims <- matrix(0,nrow=dim(ncep2.data)[1],ncol=slen)

era.file <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/grouse_mountain_ERA_800m_data.csv'
##era.file <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/snow_sites/sp_testing/spuzzum_creek_ERA_800m_data_193_121.csv'
era.data <- read.csv(era.file,header=T,as.is=T)
 
era.swe.sims <- matrix(0,nrow=dim(era.data)[1],ncol=slen)
era.snow.sims <- matrix(0,nrow=dim(era.data)[1],ncol=slen)


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
           'chilliwack_river',
           'upper_squamish',
           'tenquille_lake',
           'spuzzum_creek',
           'wahleach_lake')


for (i in seq_along(sites)) {
    site <- sites[i]
    print(site)
    ##Reanalysis 800m data
    ##ncep2.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/snow_sites/sp_testing/',site,'_NCEP2_800m_data_193_121.csv')
if (1==1) {
    ncep2.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/',site,'_NCEP2_800m_data.csv')

    ncep2.data <- read.csv(ncep2.file,header=T,as.is=T)
    ncep2.pr <- ncep2.data$Pr
    ncep2.tasmax <- ncep2.data$Tasmax
    ncep2.tasmin <- ncep2.data$Tasmin
    ncep2.tas <- ncep2.data$Tas
    ncep2.dates <- ncep2.data$Dates
}
    ##era.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/snow_sites/sp_testing/',site,'_ERA_800m_data_193_121.csv')
    era.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/',site,'_ERA_800m_data.csv')
    era.data <- read.csv(era.file,header=T,as.is=T)
    era.pr <- era.data$Pr
    era.tasmax <- era.data$Tasmax
    era.tasmin <- era.data$Tasmin
    era.tas <- era.data$Tas
    era.dates <- era.data$Dates

    coords <- get.coordinates(site)
    lat.bnds <- coords[2]
    elev <- coords[3]
    site.slope <- bc.slopes[which.min(abs(coords[1]-bc.lon)),which.min(abs(coords[2]-bc.lat))]    
    site.aspect <- bc.aspects[which.min(abs(coords[1]-bc.lon)),which.min(abs(coords[2]-bc.lat))]    


    ##coeffs <- list(a=-49.49,b=0.4128,c=2.6545,d=1.0209)
    coeffs <- list(a=-49.49,b=0.5628,c=1.5,d=1.0209)

    for (k in 1:slen) {
    print(k)
if (1==1) {
    ##ncep2.snow <- hyper.snow(ncep2.pr,ncep2.tasmax,ncep2.tasmin,coeffs)
    ncep2.results <- snow.melt(Date=ncep2.dates, precip_mm=ncep2.pr, Tmax_C=ncep2.tasmax, Tmin_C=ncep2.tasmin,
                         lat_deg=lat.bnds, slope=site.slope, aspect=site.aspect, tempHt=1, windHt=2, groundAlbedo=0.25,
                         SurfEmissiv=0.95, windSp=1, forest=0, startingSnowDepth_m=0, startingSnowDensity_kg_m3=600)
    ncep2.swe.sims[,k] <- ncep2.results$swe
    ncep2.snow.sims[,k] <- ncep2.results$snowdepth
}
    ##era.snow <- hyper.snow(era.pr,era.tasmax,era.tasmin,coeffs)
    era.results <- snow.melt(Date=era.dates, precip_mm=era.pr, Tmax_C=era.tasmax, Tmin_C=era.tasmin,
                         lat_deg=lat.bnds, slope=site.slope, aspect=site.aspect, tempHt=1, windHt=2, groundAlbedo=0.25,
                         SurfEmissiv=0.95, windSp=1, forest=0, startingSnowDepth_m=0, startingSnowDensity_kg_m3=600)
    era.swe.sims[,k] <- era.results$swe
    era.snow.sims[,k] <- era.results$snowdepth

##    Rprof(NULL)


    }                                      

    write.table(rbind(1:slen,round(ncep2.swe.sims*1000,1)),file=paste0(save.dir,site,'.ncep2.swe.',slen,'.csv'),sep=',',row.name=F,col.name=F,quote=F)
    write.table(rbind(1:slen,round(ncep2.snow.sims*100,1)),file=paste0(save.dir,site,'.ncep2.snow.',slen,'.csv'),sep=',',row.name=F,col.name=F,quote=F)

    write.table(rbind(1:slen,round(era.swe.sims*1000,1)),file=paste0(save.dir,site,'.era.swe.',slen,'.csv'),sep=',',row.name=F,col.name=F,quote=F)
    write.table(rbind(1:slen,round(era.snow.sims*100,1)),file=paste0(save.dir,site,'.era.snow.',slen,'.csv'),sep=',',row.name=F,col.name=F,quote=F)

}        
