##------------------------------------------------------------

source('/storage/home/ssobie/code/repos/winter_sports/site.coordinates.r',chdir=T)

library(ncdf4)
library(PCICt)

get_station_data <- function(ec.id,stn.dir,
                             var.name) {

  ##EC Station Data
  stn.data <- read.csv(paste(stn.dir,ec.id,'_',var.name,'.csv',sep=''),header=T,as.is=T)
  rv <- list(dates=stn.data$time,
             data=stn.data$data)  
  return(rv)
}

##-------------------------------------------------------

get_station_coordinates <- function(site) {

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

##-------------------------------------------------------

get_calibrated_parameters <- function(coords,model) {

  base.dir <- '/storage/data/projects/rci/data/winter_sports/'
  coeff.files <- list.files(path=base.dir,pattern=paste0('_hyper_snow_calibrated_parameter_',model,'_prism_TPS.nc'))
  scale.file <- coeff.files[grep('scale',coeff.files)]
  slope.file <- coeff.files[grep('slope',coeff.files)]
  freq.file <- coeff.files[grep('freq',coeff.files)]
  
  scale.nc <- nc_open(paste(base.dir,scale.file,sep=''))
  slope.nc <- nc_open(paste(base.dir,slope.file,sep=''))
  freq.nc <- nc_open(paste(base.dir,freq.file,sep=''))
  
  ##
  lon <- ncvar_get(scale.nc,'lon')
  lat <- ncvar_get(scale.nc,'lat')
  
  lon.bnds <- coords[1]
  lat.bnds <- coords[2]
  elev <- coords[3]
  
  lon.ix <- which.min(abs(lon-lon.bnds))
  lat.ix <- which.min(abs(lat-lat.bnds))
  
  scale.data <- -1 * ncvar_get(scale.nc,'scale',start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
  slope.data <- ncvar_get(slope.nc,'slope',start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
  freq.data <- ncvar_get(freq.nc,'freq',start=c(lon.ix,lat.ix,1),count=c(1,1,-1))

  nc_close(scale.nc)
  nc_close(slope.nc)
  nc_close(freq.nc)

  rv <- list(a=as.numeric(scale.data),
             b=as.numeric(slope.data),
             c=as.numeric(freq.data),
             d=1.0209)
  return(rv)  
}

##-------------------------------------------------------

get_800m_data <- function(coords,model) {

  base.dir <- paste0('/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/',model,'/')
  gcm.prism.files <- list.files(path=base.dir,pattern=paste0('gcm_prism_',model))
  pr.file <- gcm.prism.files[grep('pr_gcm',gcm.prism.files)]
  tasmax.file <- gcm.prism.files[grep('tasmax_gcm',gcm.prism.files)]
  tasmin.file <- gcm.prism.files[grep('tasmin_gcm',gcm.prism.files)]
  
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
  dates <- as.character(format(past.origin + (past.values)*86400,'%Y-%m-%d'))

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

##-------------------------------------------------------

clean <- function(input) {

  flag <- is.na(input$data)
  data <- input$data[!flag]
  dates <- input$dates[!flag]
  rv <- list(data=data,dates=dates)
  return(rv)
}

##-------------------------------------------------------

hyper_snow <- function(pr,tas,coeffs) {

        frac <- coeffs$a*(tanh(coeffs$b*(tas-coeffs$c))-coeffs$d)
        sample <- runif(length(tas),min=0,max=100)
        test <- sample > frac
        ##high.temp <- tas >= 12
        ##test[high.temp] <- TRUE
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

##-------------------------------------------------------

deterministic_snow <- function(pr,tas,coeffs) {

        frac <- coeffs$a*(tanh(coeffs$b*(tas-coeffs$c))-coeffs$d)/100

        NewSnowWatEq <- pr*frac
        R_m <- pr*(1-frac)
        rv <- list(snow=NewSnowWatEq,
                   rain=R_m)

        return(rv)
}

  
##-------------------------------------------------------

snow_vs_tas <-function(stn.id,stn.name,stn.dir,model,stn.title,y.lim,mark) {
  
  tas <- clean(get_station_data(stn.id,stn.dir,'MEAN_TEMP'))
  pr <- clean(get_station_data(stn.id,stn.dir,'ONE_DAY_PRECIPITATION'))
  snow <- clean(get_station_data(stn.id,stn.dir,'ONE_DAY_SNOW'))
  pack <- clean(get_station_data(stn.id,stn.dir,'SNOW_ON_THE_GROUND'))  

  ##Modelled snow from downscaled data
  coords <- get_station_coordinates(stn.name)

  ##Save the model series to decrease loading times
  save.dir <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/snow_sites/'
  save.file <- paste0(stn.name,'_',model,'_800m_data_series.RData')
  if (!file.exists(paste0(save.dir,save.file))) {
     data.800 <- get_800m_data(coords,model)
     save(data.800,file=paste0(save.dir,save.file))
  } else {
     load(paste0(save.dir,save.file)) 
  }

  tas.model <- data.800$tas
  pr.model <- data.800$pr
  model.dates <- data.800$dates

###Load from calibrated model parameters
  coeffs <- get_calibrated_parameters(coords,model)

  model.snow <- hyper_snow(data.800$pr,data.800$tas,coeffs)
  deter.snow <- deterministic_snow(data.800$pr,data.800$tas,coeffs)

  tas.st <- -30:29
  tas.en <- -29:40
  tas.sub <- vector(mode='list',length=70)
  for (j in 1:60) {
      tas.sub[[j]] <- tas.model >= tas.st[j] & tas.model < tas.en[j]
  }


  model.snow.matrix <- matrix(0,nrow=70,ncol=100)
  tas.matrix <- matrix(rep(tas.model,100),nrow=length(tas.model),ncol=100,byrow=F)
  for (k in 1:100) {
    tmp.snow <- hyper_snow(data.800$pr,data.800$tas,coeffs)$snow
    for (j in 1:70) {      
      tmp.count <- tmp.snow[tas.sub[[j]]]
      model.snow.matrix[j,k] <- sum(tmp.count>0.1)
    }
  }
  model.sum <- round(apply(model.snow.matrix,1,mean))
  flag <- tmp.snow > 0.1
  model.hist <- hist(tas.model[flag],breaks=-30:40,plot=F)

  print('Model')
  print(mean(tas.model[flag],na.rm=T))  
  print(sd(tas.model[flag],na.rm=T))

  snow.test.matrix <- matrix(0,nrow=length(model.snow$snow),ncol=10000)
  rain.test.matrix <- matrix(0,nrow=length(model.snow$snow),ncol=10000)
  for (k in 1:10000) {
    tmp.snow <- hyper_snow(data.800$pr,data.800$tas,coeffs)
    snow.test.matrix[,k] <- tmp.snow$snow
    rain.test.matrix[,k] <- tmp.snow$rain
  }

  snow.test.mean <- apply(snow.test.matrix,1,mean)
  rain.test.mean <- apply(rain.test.matrix,1,mean)
browser()

  print(range(tmp.snow[flag]))
  mhist <- hist(tmp.snow[flag],breaks=seq(0,150,5),plot=F)
  model.bins <- model.hist
  model.bins$counts <- round(model.sum)
  model.bins$density <- round(model.sum)/sum(model.sum,na.rm=T)

  ##Modelled snow from observed precipitation
  ##Check matching dates
  pr.match <- pr$dates %in% model.dates
  pr.data <- pr$data[pr.match]
  pr.dates <- pr$dates[pr.match]
  tas.match <- tas$dates %in% model.dates
  tas.dates <- tas$dates[tas.match]
  tas.data <- tas$data[tas.match]

  pr2.match <- pr.dates %in% tas.dates
  tas2.match <- tas.dates %in% pr.dates

  pr.data <- pr.data[pr2.match]
  tas.data <- tas.data[tas2.match]
  pr.dates <- pr.dates[pr2.match]
  tas.dates <- tas.dates[tas2.match]

  print(stn.name)


  sim.snow <- hyper_snow(pr.data,tas.data,coeffs)

  tas.sub <- vector(mode='list',length=70)
  for (j in 1:70) {
      tas.sub[[j]] <- tas.data >= tas.st[j] & tas.data < tas.en[j]
  }

  sim.snow.matrix <- matrix(0,nrow=70,ncol=100)
  tas.matrix <- matrix(rep(tas.data,100),nrow=length(tas.data),ncol=100,byrow=F)
  for (k in 1:100) {
    tmp.snow <- hyper_snow(pr.data,tas.data,coeffs)$snow
    for (j in 1:70) {      
      tmp.count <- tmp.snow[tas.sub[[j]]]
      sim.snow.matrix[j,k] <- sum(tmp.count>0.1)
    }
  }
  sim.sum <- round(apply(sim.snow.matrix,1,mean))
  flag <- tmp.snow > 0.1
  sim.hist <- hist(tas.data[flag],breaks=-30:40,plot=F)
  print('Simulation')
  print(mean(tas.data[flag],na.rm=T))  
  print(sd(tas.data[flag],na.rm=T))

  ##sim.sum <- rep(0,51)
  ##for (k in 1:70) {
  ##  tas.sub <- tas.data >= tas.st[k] & tas.data < tas.en[k]
  ##  sim.sum[k] <- mean(sim.snow$snow[tas.sub],na.rm=T)
  ##}
  ##flag <- sim.snow$snow > 0.1
  ##sim.hist <- hist(tas.data[flag],breaks=-30:40,plot=F)

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
  obs.hist <- hist(tas.obs[flag],breaks=-30:40,plot=F)
  print('Observations')
  print(mean(tas.obs[flag],na.rm=T))  
  print(sd(tas.obs[flag],na.rm=T))


  print(range(snow.data[flag]))
  snow.hist <- hist(snow.data[flag], breaks=seq(0,150,5),plot=F)
  obs.sum <- rep(0,70)
  for (k in 1:70) {
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
  plot(obs.hist,freq=T,col='lightblue',border='blue',main='Observed Snowfall',xlab='Mean Temp (\u00B0C)',ylab='Snow Frequency',
       cex=2,cex.lab=2,cex.axis=2,cex.main=2,xlim=c(-20,20),ylim=y.lim,pch='*')
  text(10,0.95*y.lim[2],stn.title,cex=2)     
  box(which='plot')
  abline(v=0)

  plot(sim.hist,freq=T,col='lightgreen',border='darkgreen',main='Simulated Snow from \nObserved Precipitation',xlab='Mean Temp (\u00B0C)',ylab='Snow Frequency',
       cex=2,cex.lab=2,cex.axis=2,cex.main=2,xlim=c(-20,20),ylim=y.lim,pch='*')
  abline(v=0)
  box(which='plot')

  plot(model.hist,freq=T,col='orange',border='darkred',main='Simulated Snow from\nModelled Precipitation',xlab='Mean Temp (\u00B0C)',ylab='Snow Frequency',
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
##  par(mar=c(5,5,5,2))
  title1 <- 'Observed Snowfall'
  title2 <- 'Simulated Snow from \nObserved Precipitation'
  title3 <-  'Simulated Snow from\nModelled Precipitation'

  plot(obs.bins,freq=F,col='lightblue',border='blue',main='',xlab='Mean Temp (\u00B0C)',ylab='Snow Frequency',
       cex=2,cex.lab=2,cex.axis=2,cex.main=2,xlim=c(-20,20),ylim=y.lim,pch='*',axes=F,yaxs='i')
 
  axis(2,at=seq(0,0.20,0.05),label=seq(0,0.20,0.05),cex.axis=2.25)
  if (mark==3) {axis(1,at=c(-10,0,10),label=c(-10,0,10),cex.axis=2.25)}

  text(-13,0.75*y.lim[2],stn.title,cex=2)     
  abline(v=0)
  if (mark==1) {
    rect(-22,0.19,22,0.24,col='white')
    text(0,0.215,title1,cex=2)
  }
  box(which='plot')  

  plot(sim.bins,freq=F,col='lightgreen',border='darkgreen',main='',xlab='Mean Temp (\u00B0C)',ylab='Snow Frequency',
       cex=2,cex.lab=2,cex.axis=2,cex.main=2,xlim=c(-20,20),ylim=y.lim,pch='*',axes=F,yaxs='i')
  abline(v=0)
  if (mark==3) {axis(1,at=c(-10,0,10),label=c(-10,0,10),cex.axis=2.25)}
  if (mark==1) {
     rect(-22,0.19,22,0.24,col='white')
     text(0,0.21,title2,cex=2)
  }
  box(which='plot')

  plot(model.bins,freq=F,col='orange',border='darkred',main='',xlab='Mean Temp (\u00B0C)',ylab='Snow Frequency',
       cex=2,cex.lab=2,cex.axis=2,cex.main=2,xlim=c(-20,20),ylim=y.lim,pch='*',axes=F,yaxs='i')
  abline(v=0)
  if (mark==1) {
     rect(-22,0.19,22,0.24,col='white')
     text(0,0.21,title3,cex=2)
  }
  if (mark==3) {axis(1,at=c(-10,0,10),label=c(-10,0,10),cex.axis=2.25)}
  box(which='plot')
  
}

if (1==0) {
  par(mar=c(5,5,5,2))
  plot(obs.bins,freq=F,col='lightblue',border='blue',main='Observed Snowfall',xlab='Mean Temp (\u00B0C)',ylab='Snow Frequency (%)',
       cex=2,cex.lab=2,cex.axis=2,cex.main=2,xlim=c(-20,20),ylim=y.lim,pch='*',axes=F)
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
stn.titles <- c('Chilliwack','Grouse','Whistler\nRoundhouse')
y.lims <- list(c(0,0.23),c(0,0.23),c(0,0.23)) ##list(c(0,90),c(0,130),c(0,110)) ##
##y.lims <- list(c(0,21),c(0,21),c(0,21)) ##list(c(0,90),c(0,130),c(0,110)) ##


base.dir <- '/storage/data/projects/rci/data/assessments/snow_model/station_data/lower_mainland/'
source('/storage/data/projects/rci/assessments/code/extract.station.daily.r',chdir=T)

model <- 'PNWNAmet'

plot.file <- paste('/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/',model,'_comparison_hyper_snow_temperature.2020.png',sep='')

###png(file=plot.file,width=6,height=6,units='in',res=600,pointsize=6,bg='white')

###par(mfrow=c(3,3))
###par(mar=c(0,0,0,0),oma=c(6,6,4,4))
###par(mgp=c(4,1.5,0))

for (i  in seq_along(stn.list)) {
  stn.id <- stn.list[i]
  stn.title <- stn.titles[i]
  stn.info <- old.get.stn.info(stn.id)
  stn.dir <- paste(base.dir,stn.info$substn,'/',sep='')
  print(stn.info$substn)
  
##  plot.file <- paste('/storage/data/projects/rci/data/winter_sports/plots/snow_temps/',model,'_',stn.info$substn,'_hyper_snow_histogram_temp.png',sep='')
##  png(file=plot.file,width=1000,height=300)
##  par(mfrow=c(1,3))

  snow_vs_tas(stn.id,stn.info$substn,stn.dir,model,stn.titles[i],y.lims[[i]],mark=i)
  ##snow.vs.tas(stn.id,stn.info$substn,stn.dir,model)

  mtext("Daily Average Temperature (\u00B0C)",side=1,outer=TRUE,cex=1.5,line=3.6)
  mtext("Density",side=2,outer=TRUE,cex=1.5,line=3.6)


##  dev.off()
}

dev.off()
}
  
