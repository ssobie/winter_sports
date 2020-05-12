##Script to compare snowdepth, SWE and snow cover between modelled values, 
##snow course sites and MODIS snow cover
library(ncdf4)
library(plotrix)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)
source('/storage/home/ssobie/code/repos/winter_sports/snow.model.functional.r',chdir=T)
source('/storage/home/ssobie/code/repos/winter_sports/site.coordinates.r',chdir=T)

hyper_snow <- function(pars,pr,tasmax,tasmin) {
           
    ##coeffs <- list(a=-48.49,b=0.7628,c=3.0,d=1.0209)
    coeffs <- pars

    tas <- (tasmax+tasmin)/2                # degrees C
    frac <- coeffs$a*(tanh(coeffs$b*(tas-coeffs$c))-coeffs$d)
    sample <- runif(length(tas),min=0,max=100)
    test <- sample > frac
    high.temp <- tas > 15
    test[high.temp] <- TRUE
    low.temp <- tas < -10
    test[low.temp] <- FALSE


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

   slen <- 100
   snow.series <- matrix(0,nrow=length(data$pr),ncol=slen)
   rain.series <- matrix(0,nrow=length(data$pr),ncol=slen)

   for (k in 1:slen) {
      snow <- hyper_snow(pars,data$pr,data$tasmax,data$tasmin)
      snow.series[,k] <- snow$swe
      rain.series[,k] <- snow$rain
   }

   snow.series.mean <- apply(snow.series,1,mean)/1000
   rain.series.mean <- apply(rain.series,1,mean)/1000
   rain.series.mean[snow.series.mean !=0] <- 0 

   snow.mean <- list(swe=snow.series.mean,rain=rain.series.mean)

   results <- snow_melt_functional(Date=data$dates, precip_mm=snow.mean, Tmax_C=data$tasmax, Tmin_C=data$tasmin,
                        lat_deg=loc[1], slope=loc[2], aspect=loc[3], tempHt=1, windHt=2, groundAlbedo=0.25, 
                        SurfEmissiv=0.95, windSp=1, forest=0, startingSnowDepth_m=0, startingSnowDensity_kg_m3=600)  

   mean(abs(data$obs[data$obs.sub] -results$swe[data$clim.sub]*1000),na.rm=T)

}

##args <- commandArgs(trailingOnly=TRUE)
##for(i in 1:length(args)){
##    eval(parse(text=args[[i]]))
##}
##model <- gcm

model <- 'PNWNAmet'
site <- 'great_bear'

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


 
##Loop over sites
##cal.sites <- c('blackwall_peak_pillow','chilliwack_river','spuzzum_creek','tenquille_lake','upper_squamish','wahleach_lake',
##               'brookmere','callaghan','dickson_lake','disappointment_lake','dog_mountain','duffey_lake','gnawed_mountain',
##               'grouse_mountain','hamilton_hill','highland_valley','klesilkwa','lightning_lake','mcgillivray_pass','nahatlatch',
##               'orchid_lake','palisade_lake','shovelnose_mountain','stave_lake','sumallo_river_west','wahleach')

##cal.sites <- c('duffey_lake','gnawed_mountain',
##               'grouse_mountain','hamilton_hill','highland_valley','klesilkwa','lightning_lake','mcgillivray_pass','nahatlatch',
##               'orchid_lake','palisade_lake','shovelnose_mountain','stave_lake','sumallo_river_west','wahleach')

##param.matrix <- matrix(0,nrow=length(cal.sites),ncol=3)


##for (k in seq_along(cal.sites)) {
##   site <- cal.sites[k]
   print(site)
   ##Reanalysis 800m data
   clim.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/',site,'_',model,'_800m_data.csv')
   clim.data <- read.csv(clim.file,header=T,as.is=T)
   clim.dates <- clim.data$Dates
   clim.tasmax <- clim.data$Tasmax
   clim.tasmin <- clim.data$Tasmin
   clim.pr <- clim.data$Pr
   
   Tmax.check <- clim.tasmax <  clim.tasmin
   clim.tasmax[Tmax.check] <- clim.tasmin[Tmax.check]+1


   ##Apply snow fall model

   coords <- get_coordinates(site)
   lat.bnds <- coords[2]
   elev <- coords[3]

   ##Snow Course Data
   course.file <- paste0('/storage/data/projects/rci/data/winter_sports/obs/snow_courses/',site,'_snow_course.csv',sep='')
   course.data <- read.csv(course.file,header=T,as.is=T)
   course.swe <- course.data[,3]
   course.dates <- format(as.Date(course.data[,1]),'%Y-%m-%d')

   ##Snow Pillow Data
##   pillow.file <- paste('/storage/data/projects/rci/data/winter_sports/obs/snow_pillow/',site,'.csv',sep='')
##   pillow.data <- read.csv(pillow.file,header=T,as.is=T)
##   flag <- is.na(pillow.data[,11])
##   pillow.dates <- format(as.Date(pillow.data[!flag,2]),'%Y-%m-%d')
##   pillow.swe <- pillow.data[!flag,11] ##mm

   obs.dates <- course.dates
   obs.swe <- course.swe

   ##--Snow Model

   site.slope <- bc.slopes[which.min(abs(coords[1]-bc.lon)),which.min(abs(coords[2]-bc.lat))]
   site.aspect <- bc.aspects[which.min(abs(coords[1]-bc.lon)),which.min(abs(coords[2]-bc.lat))]

   loc <- c(coords[2],site.slope,site.aspect)

   ##----------------------------------------------------------------------------
   ##New Calibration
   ##Split the records in half and use cross-validation
   ##Fit a,b,c sequentially

   obs.start <- max( c(min(as.Date(obs.dates)), as.Date('1980-01-01'))) 
   obs.end <- min( c(max(as.Date(obs.dates)), as.Date('2018-12-31')))
   obs.mid <- (obs.end - obs.start)/2 + obs.start
   print(paste0('Start: ',obs.start,' Mid: ',obs.mid,' End: ',obs.end))

   ##start.dates <- c(as.Date('1980-08-01'),as.Date('2009-08-01'))
   ##end.dates <- c(as.Date('2009-07-31'),as.Date('2018-07-31'))
   start.dates <- c(obs.start,obs.mid)
   end.dates <- c(obs.mid,obs.end)

   print('Check validation dates')
   
   ##Include temperature offsets for elevation biases
   tasmax.offset <- 0
   tasmin.offset <- 0

   start.pars <- list(a=-48.49,b=0.7628,c=3.0,d=1.0209)

   val.mae <- c(0,0)
   cal.pars <- matrix(NA,nrow=2,ncol=3)
   rownames(cal.pars) <- c('Val1','Val2')
   colnames(cal.pars) <- c('A','B','C')

   ##par(mfrow=c(2,3))

   for (s in seq_along(start.dates)) {

      clim.cal.sub <- clim.dates > start.dates[s] & clim.dates < end.dates[s]
      obs.cal.sub <- obs.dates > start.dates[s] & obs.dates < end.dates[s]

      clim.val.sub <- !clim.cal.sub
      obs.val.sub <- !obs.cal.sub

      clim.cal.pr <- clim.pr[clim.cal.sub]
      clim.cal.tasmax <- clim.tasmax[clim.cal.sub] + tasmax.offset
      clim.cal.tasmin <- clim.tasmin[clim.cal.sub] + tasmin.offset
      clim.cal.dates <- clim.dates[clim.cal.sub]

      obs.cal.swe <- obs.swe[obs.cal.sub]
      obs.cal.dates <- obs.dates[obs.cal.sub]

      ##Fit the paramaters to half the data    

      clim.cal.subset <- as.Date(clim.cal.dates) %in% as.Date(obs.cal.dates)
      obs.cal.subset <- as.Date(obs.cal.dates) %in% as.Date(clim.cal.dates)

      cal.data <- list(dates=clim.cal.dates,tasmax=clim.cal.tasmax,tasmin=clim.cal.tasmin,pr=clim.cal.pr,
                   obs.sub=obs.cal.subset,clim.sub=clim.cal.subset,obs=obs.cal.swe)

      ##-------------------------------------------
      clim.val.pr <- clim.pr[clim.val.sub]
      clim.val.tasmax <- clim.tasmax[clim.val.sub] + tasmax.offset
      clim.val.tasmin <- clim.tasmin[clim.val.sub] + tasmin.offset
      clim.val.dates <- clim.dates[clim.val.sub]

      obs.val.swe <- obs.swe[obs.val.sub]
      obs.val.dates <- obs.dates[obs.val.sub]

      clim.val.subset <- as.Date(clim.val.dates) %in% as.Date(obs.val.dates)
      obs.val.subset <- as.Date(obs.val.dates) %in% as.Date(clim.val.dates)

      val.data <- list(dates=clim.val.dates,tasmax=clim.val.tasmax,tasmin=clim.val.tasmin,pr=clim.val.pr,
                   obs.sub=obs.val.subset,clim.sub=clim.val.subset,obs=obs.val.swe)

      pars <- start.pars

      ##Fit A
      par.a.vals <- seq(-52,-48,0.1)
      par.a.rss <- length(par.a.vals)
      for (i in seq_along(par.a.vals)) {
         print(paste0('Parameter A increment: ',i,' of ',length(par.a.vals)))
         pars$a <- par.a.vals[i]
         par.a.rss[i] <- min.RSS(par=pars,data=cal.data,loc=loc)
      }
      pars$a <- par.a.vals[which.min(par.a.rss)]

      ##plot(par.a.vals,par.a.rss,type='l',main='Paramater A MAE')
  
      ##Fit B
      par.b.vals <- seq(0.25,1.25,0.05)
      par.b.rss <- length(par.b.vals)
      for (i in seq_along(par.b.vals)) {
         print(paste0('Parameter B increment: ',i,' of ',length(par.b.vals)))
         pars$b <- par.b.vals[i]
         par.b.rss[i] <- min.RSS(par=pars,data=cal.data,loc=loc)
      }
      pars$b <- par.b.vals[which.min(par.b.rss)]
      ##plot(par.b.vals,par.b.rss,type='l',main='Paramater B MAE')
      ##Fit C
      par.c.vals <- seq(1.5,4.5,0.05)
      par.c.rss <- length(par.c.vals)
      for (i in seq_along(par.c.vals)) {
         print(paste0('Parameter C increment: ',i,' of ',length(par.c.vals)))
         pars$c <- par.c.vals[i]
         par.c.rss[i] <- min.RSS(par=pars,data=cal.data,loc=loc)
      }
      pars$c <- par.c.vals[which.min(par.c.rss)]

      cal.pars[s,] <- c(pars$a,pars$b,pars$c)
      val.mae[s] <- min.RSS(par=pars,data=val.data,loc=loc)  

      ###plot(par.c.vals,par.c.rss,type='l',main='Paramater C MAE')
   }
   print(cal.pars)
   print('Mean')
   print(apply(cal.pars,2,mean))
   ##param.matrix[k,] <- apply(cal.pars,2,mean)
   rv <- c(site,apply(cal.pars,2,mean),1.0209)
   param.file <- paste0('/storage/data/projects/rci/data/winter_sports/calibration_parameters/',site,'_',model,'_more_limited.csv')
   write.table(rv,file=param.file,quote=FALSE,sep=',',row.name=F,col.name=F)

##}

##rv <- cbind(cal.sites,param.matrix,rep(1.0209,length(cal.sites)))

