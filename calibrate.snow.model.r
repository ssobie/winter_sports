##Script to compare snowdepth, SWE and snow cover between modelled values, 
##snow course sites and MODIS snow cover
library(ncdf4)
library(plotrix)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)
source('/storage/home/ssobie/code/repos/winter_sports/snow.model.functional.r',chdir=T)


hyper_snow <- function(pars,pr,tasmax,tasmin) {
           
    ##coeffs <- list(a=-48.49,b=0.7628,c=3.0,d=1.0209)
    coeffs <- pars

    tas <- (tasmax+tasmin)/2                # degrees C
    frac <- coeffs$a*(tanh(coeffs$b*(tas-coeffs$c))-coeffs$d)
    sample <- runif(length(tas),min=0,max=100)
    test <- sample > frac
    snow.type <- rep(TRUE,length(tas))
    snow.type[test] <- FALSE

    NewSnowWatEq <- pr
    NewSnowWatEq[!snow.type] <- 0
    R_m <- pr
    R_m[snow.type] <- 0
    rv <- list(swe=NewSnowWatEq,
               rain=R_m)
    return(rv)
}

min.RSS <- function(par,data,loc) {

   slen <- 10
   snow.series <- matrix(0,nrow=length(data$pr),ncol=slen)
   rain.series <- matrix(0,nrow=length(data$pr),ncol=slen)
   for (k in 1:slen) {
      snow <- hyper_snow(pars,data$pr,data$tasmax,data$tasmin)
      snow.series[,k] <- snow$swe
      rain.series[,k] <- snow$rain
   }

   snow.series.mean <- apply(snow.series,1,mean)
   rain.series.mean <- apply(snow.series,1,mean)
   rain.series.mean[snow.series.mean !=0] <- 0 

   snow.mean <- list(swe=snow.series.mean*0.001,rain=rain.series.mean*0.001)

   results <- snow.melt(Date=data$dates, precip_mm=snow.mean, Tmax_C=data$tasmax, Tmin_C=data$tasmin,
                        lat_deg=loc[1], slope=loc[2], aspect=loc[3], tempHt=1, windHt=2, groundAlbedo=0.25,                                   SurfEmissiv=0.95, windSp=1, forest=0, startingSnowDepth_m=0, startingSnowDensity_kg_m3=600)  

   mean(abs(data$obs[data$obs.sub] -results$swe[data$clim.sub]*1000),na.rm=T)

}



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

##North Shore Sites
site <- 'spuzzum_creek'
site.name <- 'Spuzzum Creek'

model <- 'ERA'

##Loop over sites
model.dir <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sims/'

print(site)
##Reanalysis 800m data
clim.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/',site,'_',model,'_800m_data.csv')
clim.data <- read.csv(clim.file,header=T,as.is=T)
clim.dates <- clim.data$Dates

##Apply snow fall model

coords <- get.coordinates(site)
lat.bnds <- coords[2]
elev <- coords[3]

##Snow Course Data
##course.file <- paste0('/storage/data/projects/rci/data/winter_sports/obs/snow_courses/',site,'_snow_course.csv',sep='')
##course.data <- read.csv(course.file,header=T,as.is=T)
##course.dates <- format(as.Date(course.data[,1]),'%Y-%m-%d')

##Snow Pillow Data
pillow.file <- paste('/storage/data/projects/rci/data/winter_sports/obs/snow_pillow/',site,'_asp.csv',sep='')
pillow.data <- read.csv(pillow.file,header=T,as.is=T)
flag <- is.na(pillow.data[,11])
pillow.dates <- format(as.Date(pillow.data[!flag,2]),'%Y-%m-%d')
pillow.swe <- pillow.data[!flag,11] ##mm

obs.dates <- pillow.dates
obs.swe <- pillow.swe

swe.sims <- read.csv(paste0(model.dir,site,'.',tolower(model),'.swe.1001.csv'),header=T,as.is=T)
swe.mean <- apply(swe.sims,1,mean,na.rm=T)
rm(swe.sims)

##--Snow Model

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

site.slope <- bc.slopes[which.min(abs(coords[1]-bc.lon)),which.min(abs(coords[2]-bc.lat))]
site.aspect <- bc.aspects[which.min(abs(coords[1]-bc.lon)),which.min(abs(coords[2]-bc.lat))]

loc <- c(coords[2],site.slope,site.aspect)

##--------------------------------------------------------------------------------
##Existing Calibration Method
if(1==0) {

clim.subset <- as.Date(clim.dates) %in% as.Date(course.dates)
course.subset <- as.Date(course.dates) %in% as.Date(clim.dates)

##start.pars <- list(a=-48.49,b=0.7628,c=3.0,d=1.0209)
pars <-  list(a=-43.49,b=0.7628,c=0.0,d=1.0209)
data <- list(dates=clim.dates,tasmax=clim.tasmax,tasmin=clim.tasmin,pr=clim.pr,
             course.sub=course.subset,clim.sub=clim.subset,course=course.swe)
loc <- c(coords[2],site.slope,site.aspect)

par.tmp <- seq(-50,-40,0.05)
##par.tmp <- seq(0.0,2.0,0.05)
##par.tmp <- seq(2,6.0,0.05)

par.tmp.rss <- length(par.tmp)
for (i in seq_along(par.tmp)) {
  print(paste0('Parameter increment: ',i,' of ',length(par.tmp)))
  pars$a <- par.tmp[i]
  par.tmp.rss[i] <- min.RSS(par=pars,data=data,loc=loc)
}

plot(par.tmp,par.tmp.rss,type='l')

##optim.result <- optim(par=pars,fn=min.RSS,data=data,loc=loc)
##                      ##upper=upper.par,lower=lower.par,method="L-BFGS-B")
##lines(as.Date(clim.dates),results$swe*1000,col='red',lwd=2)
}
##----------------------------------------------------------------------------
##New Calibration
##Split the records in half and use cross-validation
##Fit a,b,c sequentially

start.dates <- c(as.Date('1980-08-01'),as.Date('2009-08-01'))
end.dates <- c(as.Date('2009-07-31'),as.Date('2018-07-31'))
print('Check validation dates')

##Include temperature offsets for elevation biases
tasmax.offset <- 0
tasmin.offset <- 0

start.pars <- list(a=-48.49,b=0.7628,c=3.0,d=1.0209)

val.mae <- c(0,0)
cal.pars <- matrix(NA,nrow=2,ncol=3)
rownames(cal.pars) <- c('Val1','Val2')
colnames(cal.pars) <- c('A','B','C')

par(mfrow=c(2,3))

for (s in seq_along(start.dates)) {

  clim.cal.sub <- clim.dates > start.dates[s] & clim.dates < end.dates[s]
  obs.cal.sub <- obs.dates > start.dates[s] & obs.dates < end.dates[s]

  clim.val.sub <- !clim.cal.sub
  obs.val.sub <- !obs.cal.sub

  clim.cal.pr <- clim.data$Pr[clim.cal.sub]
  clim.cal.tasmax <- clim.data$Tasmax[clim.cal.sub] + tasmax.offset
  clim.cal.tasmin <- clim.data$Tasmin[clim.cal.sub] + tasmin.offset
  clim.cal.dates <- clim.dates[clim.cal.sub]

  obs.cal.swe <- obs.swe[obs.cal.sub]
  obs.cal.dates <- obs.dates[obs.cal.sub]

  ##Fit the paramaters to half the data    

  clim.cal.subset <- as.Date(clim.cal.dates) %in% as.Date(obs.cal.dates)
  obs.cal.subset <- as.Date(obs.cal.dates) %in% as.Date(clim.cal.dates)

  cal.data <- list(dates=clim.cal.dates,tasmax=clim.cal.tasmax,tasmin=clim.cal.tasmin,pr=clim.cal.pr,
                   obs.sub=obs.cal.subset,clim.sub=clim.cal.subset,obs=obs.cal.swe)

  ##-------------------------------------------
  clim.val.pr <- clim.data$Pr[clim.val.sub]
  clim.val.tasmax <- clim.data$Tasmax[clim.val.sub] + tasmax.offset
  clim.val.tasmin <- clim.data$Tasmin[clim.val.sub] + tasmin.offset
  clim.val.dates <- clim.dates[clim.val.sub]

  obs.val.swe <- obs.swe[obs.val.sub]
  obs.val.dates <- obs.dates[obs.val.sub]

  clim.val.subset <- as.Date(clim.val.dates) %in% as.Date(obs.val.dates)
  obs.val.subset <- as.Date(obs.val.dates) %in% as.Date(clim.val.dates)

  val.data <- list(dates=clim.val.dates,tasmax=clim.val.tasmax,tasmin=clim.val.tasmin,pr=clim.val.pr,
                   obs.sub=obs.val.subset,clim.sub=clim.val.subset,obs=obs.val.swe)

  pars <- start.pars

  ##Fit A
  par.a.vals <- seq(-52,-40,0.1)
  par.a.rss <- length(par.a.vals)
  for (i in seq_along(par.a.vals)) {
    print(paste0('Parameter A increment: ',i,' of ',length(par.a.vals)))
    pars$a <- par.a.vals[i]
    par.a.rss[i] <- min.RSS(par=pars,data=cal.data,loc=loc)
  }
  pars$a <- par.a.vals[which.min(par.a.rss)]

  plot(par.a.vals,par.a.rss,type='l',main='Paramater A MAE')
  
  ##Fit B
  par.b.vals <- seq(0.05,1.0,0.05)
  par.b.rss <- length(par.b.vals)
  for (i in seq_along(par.b.vals)) {
    print(paste0('Parameter B increment: ',i,' of ',length(par.b.vals)))
    pars$b <- par.b.vals[i]
    par.b.rss[i] <- min.RSS(par=pars,data=cal.data,loc=loc)
  }
  pars$b <- par.b.vals[which.min(par.b.rss)]
  plot(par.b.vals,par.b.rss,type='l',main='Paramater B MAE')
  ##Fit C
  par.c.vals <- seq(0.0,4.5,0.05)
  par.c.rss <- length(par.c.vals)
  for (i in seq_along(par.c.vals)) {
    print(paste0('Parameter C increment: ',i,' of ',length(par.c.vals)))
    pars$c <- par.c.vals[i]
    par.c.rss[i] <- min.RSS(par=pars,data=cal.data,loc=loc)
  }
  pars$c <- par.c.vals[which.min(par.c.rss)]

  cal.pars[s,] <- c(pars$a,pars$b,pars$c)
  val.mae[s] <- min.RSS(par=pars,data=val.data,loc=loc)  

  plot(par.c.vals,par.c.rss,type='l',main='Paramater C MAE')


}
print(cal.pars)
print('Mean')
print(apply(cal.pars,2,mean))
