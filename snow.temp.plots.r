##------------------------------------------------------------


library(ncdf4)
library(PCICt)

get.station.data <- function(ec.id,stn.dir,
                             var.name) {

  ##EC Station Data
  stn.data <- read.csv(paste(stn.dir,ec.id,'_',var.name,'.csv',sep=''),header=T,as.is=T)
  rv <- list(dates=stn.data$time,
             data=stn.data$data)  
  return(rv)
}

get.station.coordinates <- function(site) {

  coordinates <- list(abbotsford=c(-122.360,49.02528,59),
                      agassiz=c(-121.7597,49.2425,15),
                      alouette_lake=c(-122.48333,49.28333,117),
                      alta_lake=c(-122.95,50.15,50),
                      burnaby_simon_fraser_u=c(-122.918,49.27839,365),
                      chilliwack=c(-121.9247,49.17216,11),
                      coquitlam_lake=c(-122.800,49.36667,161),
                      grouse_mountain=c(-123.0806,49.36806,1128),
                      hope_slide=c(-121.2364,49.275,687.8),
                      mission_west_abbey=c(-122.2708,49.1525,49),
                      pemberton_airport=c(-122.7378,50.3025,204),
                      stave_falls=c(-121.3667,49.233,110),
                      whistler_roundhouse=c(-122.95,50.0667,1835),
                      whistler=c(-122.9548,50.12889,658))                      
  rv <- coordinates[[as.character(site)]]
  return(rv)
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
                      hamilton_hill=c(-120.7955805,49.4988027,1477))
  
  rv <- coordinates[[site]]
  return(rv)
}


get.800m.data <- function(coords,model) {

  
  base.dir <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/',model,'/')
  
  pr.file <- paste0('pr_gcm_prism_',model,'_1979-2016.nc')
  tasmax.file <- paste0('tasmax_gcm_prism_',model,'_1979-2016.nc')
  tasmin.file <- paste0('tasmin_gcm_prism_',model,'_1979-2016.nc')
  
  pr.nc <- nc_open(paste(base.dir,pr.file,sep=''))
  tasmax.nc <- nc_open(paste(base.dir,tasmax.file,sep=''))
  tasmin.nc <- nc_open(paste(base.dir,tasmin.file,sep=''))
  
  time.atts <- ncatt_get(pr.nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  
  time.start <- as.Date(strsplit(time.units, ' ')[[1]][3])
  past.origin <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                          cal=time.calendar)
  
  past.values <- ncvar_get(pr.nc,'time')
  dates <- as.character(format(past.origin + (past.values)*86400/24,'%Y-%m-%d'))

  ##
  lon <- ncvar_get(pr.nc,'lon')
  lat <- ncvar_get(pr.nc,'lat')
  
  lon.bnds <- coords[1]
  lat.bnds <- coords[2]
  elev <- coords[3]
  
  lon.ix <- which.min(abs(lon-lon.bnds))
  lat.ix <- which.min(abs(lat-lat.bnds))
  
  print(lon[lon.ix])
  print(lat[lat.ix])

  pr.data <- ncvar_get(pr.nc,'pr',start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
  tasmax.data <- ncvar_get(tasmax.nc,'tasmax',start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
  tasmin.data <- ncvar_get(tasmin.nc,'tasmin',start=c(lon.ix,lat.ix,1),count=c(1,1,-1))


  nc_close(pr.nc)
  nc_close(tasmax.nc)
  nc_close(tasmin.nc)

  
  flag <- which(tasmax.data<tasmin.data)
  tas.diff <- tasmax.data-tasmin.data
  tasmax.data[flag] <- tasmax.data[flag]+max(tas.diff)+0.1
  
  tas.data <- (tasmax.data+tasmin.data)/2

  rv <- list(pr=pr.data,
             tasmax=tasmax.data,
             tasmin=tasmin.data,
             tas=tas.data,
             dates=dates)
  return(rv)  
}

clean <- function(input) {

  flag <- is.na(input$data)
  data <- input$data[!flag]
  dates <- input$dates[!flag]
  rv <- list(data=data,dates=dates)
  return(rv)
}



hyper.snow <- function(pr,tas,coeffs) {

        frac <- coeffs$a*(tanh(coeffs$b*(tas-coeffs$c))-coeffs$d)
        sample <- runif(length(tas),min=0,max=100)
        test <- sample > frac
        high.temp <- tas > 12
        test[high.temp] <- TRUE
        snow.type <- rep(TRUE,length(tas))
        snow.type[test] <- FALSE

        NewSnowWatEq <- pr
        NewSnowWatEq[!snow.type] <- 0
        R_m <- pr
        R_m[snow.type] <- 0
        rv <- list(snow=NewSnowWatEq,
                   rain=R_m)
        return(rv)
}
  

snow.vs.tas <-function(stn.id,stn.name,stn.dir,model,stn.title,y.lim) {
  
  tas <- clean(get.station.data(stn.id,stn.dir,'MEAN_TEMP'))
  pr <- clean(get.station.data(stn.id,stn.dir,'ONE_DAY_PRECIPITATION'))
  snow <- clean(get.station.data(stn.id,stn.dir,'ONE_DAY_SNOW'))
  pack <- clean(get.station.data(stn.id,stn.dir,'SNOW_ON_THE_GROUND'))  

  ##Modelled snow from downscaled data
  coords <- get.station.coordinates(stn.name)

  data.800 <- get.800m.data(coords,model)
  tas.model <- data.800$tas
  pr.model <- data.800$pr

  ##coeffs <- list(a=-49.49,b=0.4128,c=2.6545,d=1.0209)
  coeffs <- list(a=-49.49,b=0.5628,c=1.5,d=1.0209)
  model.snow <- hyper.snow(data.800$pr,data.800$tas,coeffs)
   
  ##Modelled snow from observed precipitation
  ##Check matching dates
  tas.match <- tas$dates %in% pr$dates
  tas.sim <- tas$data[tas.match]
  tas.dates <- tas$dates[tas.match]

  pr.match <- pr$dates %in% tas$dates
  pr.sim <- pr$data[pr.match]
  pr.dates <- pr$dates[pr.match]

  sim.snow <- hyper.snow(pr.sim,tas.sim,coeffs)

  ##Observed Snow
  snow.match <- snow$dates %in% tas$dates
  snow.data <- snow$data[snow.match]
  snow.dates <- snow$dates[snow.match]

  tas.match <- tas$dates %in% snow$dates
  tas.data <- tas$data[tas.match]
  tas.dates <- tas$dates[tas.match]

  ##y.lim <- range(c(range(snow.data),range(sim.snow),range(model.snow)))

  par(mar=c(5,5,5,2))
  plot(tas.data,snow.data,col='blue',main='Observed Snowfall',xlab='Mean Temp (\u00B0C)',ylab='Snow (cm)',
       cex=2,cex.lab=2,cex.axis=2,cex.main=2,xlim=c(-20,20),ylim=y.lim,pch='*')
  text(10,0.95*y.lim[2],stn.title,cex=2)     
  box(which='plot')
  abline(v=0)

  plot(tas.sim,sim.snow$snow,col='green',main='Simulated Snow from \nObserved Precipitation',xlab='Mean Temp (\u00B0C)',ylab='Snow (cm)',
       cex=2,cex.lab=2,cex.axis=2,cex.main=2,xlim=c(-20,20),ylim=y.lim,pch='*')
  abline(v=0)
  box(which='plot')

  plot(tas.model,model.snow$snow,col='red',main='Simulated Snow from\nModelled Precipitation',xlab='Mean Temp (\u00B0C)',ylab='Snow (cm)',
       cex=2,cex.lab=2,cex.axis=2,cex.main=2,xlim=c(-20,20),ylim=y.lim,pch='*')
  abline(v=0)
  box(which='plot')

}



##************************************************************************************

##ANUSPLIN Test
if (1==1) {
stn.list <- c('1101530','1100120','1107680','1100030',
              '1100360','1101890','1040390','1105192','1101158',
              '1108906','1113581','1048898','1086082','1105658')

stn.list <- c('1101530','1105658','1108906')
stn.titles <- c('Chilliwack','Grouse','Whistler')
y.lims <- list(c(0,100),c(0,150),c(0,120))


base.dir <- '/storage/data/projects/rci/data/assessments/snow_model/station_data/lower_mainland/'
source('/storage/data/projects/rci/assessments/code/extract.station.daily.r',chdir=T)

model <- 'ERA'

plot.file <- paste('/storage/data/projects/rci/data/winter_sports/plots/snow_temps/',model,'_comparison_hyper_snow_temperature2.png',sep='')
png(file=plot.file,width=1000,height=1000)
par(mfrow=c(3,3))

for (i  in seq_along(stn.list)) {
  stn.id <- stn.list[i]
  stn.title <- stn.titles[i]
  stn.info <- old.get.stn.info(stn.id)
  stn.dir <- paste(base.dir,stn.info$substn,'/',sep='')
  print(stn.info$substn)
  
  ##plot.file <- paste('/storage/data/projects/rci/data/winter_sports/plots/snow_temps/',model,'_',stn.info$substn,'_hyper_snow_temperature2.png',sep='')
  ##png(file=plot.file,width=1000,height=300)
  ##par(mfrow=c(1,3))

  snow.vs.tas(stn.id,stn.info$substn,stn.dir,model,stn.title,y.lims[[i]])
  ##snow.vs.tas(stn.id,stn.info$substn,stn.dir,model)

##  dev.off()
}

dev.off()
}
  
