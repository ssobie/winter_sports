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
  model.dates <- data.800$dates

  ##coeffs <- list(a=-49.49,b=0.4128,c=2.6545,d=1.0209)
  coeffs <- list(a=-49.49,b=0.5628,c=1.5,d=1.0209)
  model.snow <- hyper.snow(data.800$pr,data.800$tas,coeffs)

  tas.st <- -30:19
  tas.en <- -29:20
  tas.sub <- vector(mode='list',length=50)
  for (j in 1:50) {
      tas.sub[[j]] <- tas.model >= tas.st[j] & tas.model < tas.en[j]
  }

  model.snow.matrix <- matrix(0,nrow=50,ncol=100)
  tas.matrix <- matrix(rep(tas.model,100),nrow=length(tas.model),ncol=100,byrow=F)
  for (k in 1:100) {
    tmp.snow <- hyper.snow(data.800$pr,data.800$tas,coeffs)$snow
    for (j in 1:50) {      
      tmp.count <- tmp.snow[tas.sub[[j]]]
      model.snow.matrix[j,k] <- sum(tmp.count>0.1)
    }
  }
  model.sum <- round(apply(model.snow.matrix,1,mean))
  flag <- tmp.snow > 0.1
  model.hist <- hist(tas.model[flag],breaks=-30:20,plot=F)

  print('Model')
  print(mean(tas.model[flag],na.rm=T))  
  print(sd(tas.model[flag],na.rm=T))

  print(range(tmp.snow[flag]))
  mhist <- hist(tmp.snow[flag],breaks=seq(0,150,5),plot=F)
  model.bins <- model.hist
  model.bins$counts <- round(model.sum)
  model.bins$density <- round(model.sum)/sum(model.sum,na.rm=T)

  ##Modelled snow from observed precipitation
  ##Check matching dates
  pr.match <- pr$dates %in% model.dates
  pr.dates <- pr$dates[pr.match]
  tas.match <- tas$dates %in% model.dates
  tas.dates <- tas$dates[tas.match]

  pr2.match <- pr.match[pr.dates %in% tas.dates]
  tas2.match <- tas.match[tas.dates %in% pr.dates]

  pr.data <- pr$data[pr2.match]
  tas.data <- tas$data[tas2.match]
  pr.dates <- pr$dates[pr2.match]
  tas.dates <- tas$dates[tas2.match]

  sim.snow <- hyper.snow(pr.data,tas.data,coeffs)

  tas.sub <- vector(mode='list',length=50)
  for (j in 1:50) {
      tas.sub[[j]] <- tas.data >= tas.st[j] & tas.data < tas.en[j]
  }

  sim.snow.matrix <- matrix(0,nrow=50,ncol=100)
  tas.matrix <- matrix(rep(tas.data,100),nrow=length(tas.data),ncol=100,byrow=F)
  for (k in 1:100) {
    tmp.snow <- hyper.snow(pr.data,tas.data,coeffs)$snow
    for (j in 1:50) {      
      tmp.count <- tmp.snow[tas.sub[[j]]]
      sim.snow.matrix[j,k] <- sum(tmp.count>0.1)
    }
  }
  sim.sum <- round(apply(sim.snow.matrix,1,mean))
  flag <- tmp.snow > 0.1
  sim.hist <- hist(tas.data[flag],breaks=-30:20,plot=F)
  print('Simulation')
  print(mean(tas.data[flag],na.rm=T))  
  print(sd(tas.data[flag],na.rm=T))

  ##sim.sum <- rep(0,51)
  ##for (k in 1:50) {
  ##  tas.sub <- tas.data >= tas.st[k] & tas.data < tas.en[k]
  ##  sim.sum[k] <- mean(sim.snow$snow[tas.sub],na.rm=T)
  ##}
  ##flag <- sim.snow$snow > 0.1
  ##sim.hist <- hist(tas.data[flag],breaks=-30:20,plot=F)

  sim.bins <- sim.hist
  sim.bins$counts <- round(sim.sum)
  sim.bins$density <- round(sim.sum)/sum(sim.sum,na.rm=T)

  ##yr.fac <- as.factor(format(as.Date(pr.dates),'%Y'))
  ##mn.fac <- as.factor(format(as.Date(pr.dates),'%m'))
  ##mn.sum <- tapply(sim.snow$snow,list(yr.fac,mn.fac),sum,na.rm=T)
  ##sim.mn <- mn.sum[,1]
  ##sim.avg <- round(apply(mn.sum,2,mean,na.rm=T),1)


  ##yr.fac <- as.factor(format(as.Date(model.dates),'%Y'))
  ##mn.fac <- as.factor(format(as.Date(model.dates),'%m'))
  ##mn.sum <- tapply(model.snow$snow,list(yr.fac,mn.fac),sum,na.rm=T)
  ##model.avg <- round(apply(mn.sum,2,mean,na.rm=T),1)
  ##model.mn <- mn.sum[,1]

  ##----------------------------------------------------
  ##Observed Snow 
  snow.match <- snow$dates %in% tas.dates
  tas.match <- tas$dates %in% snow$dates[snow.match]
  tas.obs <- tas$data[tas.match]
  snow.data <- snow$data[snow.match]
  snow.dates <- snow$dates[snow.match]

  flag <- snow.data > 0.1
  obs.hist <- hist(tas.obs[flag],breaks=-30:20,plot=F)
  print('Observations')
  print(mean(tas.obs[flag],na.rm=T))  
  print(sd(tas.obs[flag],na.rm=T))


  print(range(snow.data[flag]))
  snow.hist <- hist(snow.data[flag], breaks=seq(0,150,5),plot=F)
  obs.sum <- rep(0,50)
  for (k in 1:50) {
    tas.sub <- tas.obs >= tas.st[k] & tas.obs < tas.en[k]
    obs.sum[k] <- sum(snow.data[tas.sub] > 0.1)
  }
  obs.sum[is.na(obs.sum)] <- 0
  obs.bins <- obs.hist
  obs.bins$counts <- round(obs.sum)
  obs.bins$density <- round(obs.sum)/sum(obs.sum,na.rm=T)

  ##yr.fac <- as.factor(format(as.Date(snow.dates),'%Y'))
  ##mn.fac <- as.factor(format(as.Date(snow.dates),'%m'))
  ##mn.sum <- tapply(snow.data,list(yr.fac,mn.fac),sum,na.rm=T)
  ##obs.avg <- round(apply(mn.sum,2,mean,na.rm=T),1)
  ##obs.mn <- mn.sum[,1]
  ##print(rbind(obs.avg,sim.avg,model.avg))



##  plot(names(obs.mn),obs.mn,xlim=c(1978,2016),ylim=c(0,200),type='l',lwd=4,col='blue')
##  lines(names(sim.mn),sim.mn,lwd=4,col='green')
##  lines(names(model.mn),model.mn,lwd=4,col='red')

  ##y.lim <- range(c(range(snow.data),range(sim.snow),range(model.snow)))

if (1==0) {
  par(mar=c(5,5,5,2))
  plot(obs.hist,freq=F,col='lightblue',border='blue',main='Observed Snowfall',xlab='Mean Temp (\u00B0C)',ylab='Snow Frequency',
       cex=2,cex.lab=2,cex.axis=2,cex.main=2,xlim=c(-20,20),ylim=y.lim,pch='*')
  text(10,0.95*y.lim[2],stn.title,cex=2)     
  box(which='plot')
  abline(v=0)

  plot(sim.hist,freq=F,col='lightgreen',border='darkgreen',main='Simulated Snow from \nObserved Precipitation',xlab='Mean Temp (\u00B0C)',ylab='Snow Frequency',
       cex=2,cex.lab=2,cex.axis=2,cex.main=2,xlim=c(-20,20),ylim=y.lim,pch='*')
  abline(v=0)
  box(which='plot')

  plot(model.hist,freq=F,col='orange',border='darkred',main='Simulated Snow from\nModelled Precipitation',xlab='Mean Temp (\u00B0C)',ylab='Snow Frequency',
       cex=2,cex.lab=2,cex.axis=2,cex.main=2,xlim=c(-20,20),ylim=y.lim,pch='*')
  abline(v=0)
  box(which='plot')
}

if (1==0) {
  par(mar=c(5,5,5,2))
  flag <- snow.data > 0.1
  plot(tas.obs[flag],snow.data[flag],col='blue',main='Observed Snowfall',xlab='Mean Temp (\u00B0C)',ylab='Snow (cm)',
       cex=2.2,cex.lab=2,cex.axis=2,cex.main=2,xlim=c(-20,20),ylim=y.lim,pch='*')
  text(10,0.95*y.lim[2],stn.title,cex=2)     
  box(which='plot')
  abline(v=0)

  flag <- sim.snow$snow > 0.1
  plot(tas.data[flag],sim.snow$snow[flag],col='green',main='Simulated Snow from \nObserved Precipitation',xlab='Mean Temp (\u00B0C)',ylab='Snow (cm)',
       cex=2.2,cex.lab=2,cex.axis=2,cex.main=2,xlim=c(-20,20),ylim=y.lim,pch='*')
  abline(v=0)
  box(which='plot')

  flag <- model.snow$snow > 0.1
  plot(tas.model[flag],model.snow$snow[flag],col='red',main='Simulated Snow from\nModelled Precipitation',xlab='Mean Temp (\u00B0C)',ylab='Snow (cm)',
       cex=2.2,cex.lab=2,cex.axis=2,cex.main=2,xlim=c(-20,20),ylim=y.lim,pch='*')
  abline(v=0)
  box(which='plot')
}


if (1==1) {
  par(mar=c(5,5,5,2))
  plot(obs.bins,freq=F,col='lightblue',border='blue',main='Observed Snowfall',xlab='Mean Temp (\u00B0C)',ylab='Snow Frequency (%)',
       cex=2,cex.lab=2,cex.axis=2,cex.main=2,xlim=c(-20,20),ylim=y.lim,pch='*')
  text(10,0.95*y.lim[2],stn.title,cex=2)     
  box(which='plot')
  abline(v=0)

  plot(sim.bins,freq=F,col='lightgreen',border='darkgreen',main='Simulated Snow from \nObserved Precipitation',xlab='Mean Temp (\u00B0C)',ylab='Snow Frequency (%)',
       cex=2,cex.lab=2,cex.axis=2,cex.main=2,xlim=c(-20,20),ylim=y.lim,pch='*')
  abline(v=0)
  box(which='plot')

  plot(model.bins,freq=F,col='orange',border='darkred',main='Simulated Snow from\nModelled Precipitation',xlab='Mean Temp (\u00B0C)',ylab='Snow Frequency (%)',
       cex=2,cex.lab=2,cex.axis=2,cex.main=2,xlim=c(-20,20),ylim=y.lim,pch='*')
  abline(v=0)
  box(which='plot')
}

if (1==0) {
  par(mar=c(5,5,5,2))
  plot(obs.bins,freq=F,col='lightblue',border='blue',main='Observed Snowfall',xlab='Mean Temp (\u00B0C)',ylab='Snow Frequency (%)',
       cex=2,cex.lab=2,cex.axis=2,cex.main=2,xlim=c(-20,20),ylim=y.lim,pch='*')
  text(10,0.95*y.lim[2],stn.title,cex=2)     
  box(which='plot')
  abline(v=0)

  plot(model.bins,freq=F,col='orange',border='darkred',main='Simulated Snow from\nModelled Precipitation',xlab='Mean Temp (\u00B0C)',ylab='Snow Frequency (%)',
       cex=2,cex.lab=2,cex.axis=2,cex.main=2,xlim=c(-20,20),ylim=y.lim,pch='*')
  abline(v=0)
  box(which='plot')

  plot(snow.hist,freq=F,col='lightgreen',border='darkgreen',main='Simulated Snow from \nObserved Precipitation',xlab='Mean Temp (\u00B0C)',ylab='Snow Frequency (%)',
       cex=2,cex.lab=2,cex.axis=2,cex.main=2,xlim=c(0,70),ylim=c(0,0.25))
  abline(v=0)
  box(which='plot')


}


print('Sim-Obs')
print(sum(((sim.bins$counts+1)-(obs.bins$counts+1))^2 / (obs.bins$counts+1)))
print('Model-Obs')
print(sum(((model.bins$counts+1)-(obs.bins$counts+1))^2 / (obs.bins$counts+1)))

print(chisq.test(cbind(obs.bins$counts+1,model.bins$counts+1)))



}



##************************************************************************************

##ANUSPLIN Test
if (1==1) {
stn.list <- c('1101530','1100120','1107680','1100030',
              '1100360','1101890','1040390','1105192','1101158',
              '1108906','1113581','1048898','1086082','1105658')

stn.list <- c('1101530','1105658','1108906')
stn.titles <- c('Chilliwack','Grouse','Whistler')
y.lims <- list(c(0,0.21),c(0,0.21),c(0,0.21)) ##list(c(0,90),c(0,130),c(0,110)) ##


base.dir <- '/storage/data/projects/rci/data/assessments/snow_model/station_data/lower_mainland/'
source('/storage/data/projects/rci/assessments/code/extract.station.daily.r',chdir=T)

model <- 'ERA'

plot.file <- paste('/storage/data/projects/rci/data/winter_sports/plots/snow_temps/',model,'_comparison_hyper_snow_temperature5.png',sep='')
png(file=plot.file,width=1000,height=1000)
par(mfrow=c(3,3))

for (i  in seq_along(stn.list)) {
  stn.id <- stn.list[i]
  stn.title <- stn.titles[i]
  stn.info <- old.get.stn.info(stn.id)
  stn.dir <- paste(base.dir,stn.info$substn,'/',sep='')
  print(stn.info$substn)
  
##  plot.file <- paste('/storage/data/projects/rci/data/winter_sports/plots/snow_temps/',model,'_',stn.info$substn,'_hyper_snow_histogram_temp.png',sep='')
##  png(file=plot.file,width=1000,height=300)
##  par(mfrow=c(1,3))

  snow.vs.tas(stn.id,stn.info$substn,stn.dir,model,stn.titles[i],y.lims[[i]])
  ##snow.vs.tas(stn.id,stn.info$substn,stn.dir,model)

##  dev.off()
}

dev.off()
}
  
