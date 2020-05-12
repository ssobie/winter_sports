##Script to plot time series of frost free days for Vancouver Intl.

library(ncdf4)
library(PCICt)
library(raster)
library(rgdal)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)


##--------------------------------------------------------------------------
##Start, Peak, End Dates, and Lengths

get_snow_season_dates <- function(input.snow,dates) {

  years <- format(as.Date(dates),'%Y')
  jdays <- format(as.Date(dates),'%j')
  jdays.list <- tapply(jdays,as.factor(years),list)
  peaks <- tapply(input.snow,as.factor(years),function(x){which.max(x[1:300])})
  lows <- tapply(input.snow,as.factor(years),function(x){which.min(x[40:300])})

  sim.snow <- input.snow
  sim.snow[sim.snow>0] <- 100

  patterns <- rle(sim.snow)
  pt     <- patterns$lengths
  pt.sum <- cumsum(patterns$lengths)
  all.zeroes <- which(patterns$values==100)
  zeroes <- all.zeroes[patterns$lengths[all.zeroes] > 100]
  top.ix <- pt.sum[zeroes]
  bot.ix <- pt[zeroes]

  pt.diff <- c(pt[1],top.ix[1], top.ix[-1] - bot.ix[-1], tail(top.ix,1))
  ix <- sort(unique(c(pt.diff,top.ix)))
  if (ix[1] > 100) {
     ix <- sort(unique(c(1,pt.diff,top.ix)))
  }
  ##plot(as.Date(dates),input.snow,xlim=c(as.Date('1967-01-01'),as.Date('1980-01-01')))
  ##abline(v=as.Date(dates[ix]),col='red')
  len <- length(ix)

  lengths <- ix[seq(2,len,2)] - ix[seq(1,len-1,2)]
  ##starts <- as.numeric(format(as.Date(dates[ix[seq(1,len-1,2)]]),'%j'))
  ##ends <- as.numeric(format(as.Date(dates[ix[seq(2,len,2)]]),'%j'))
  starts <- as.Date(dates[ix[seq(1,len-1,2)]])
  ends <- as.Date(dates[ix[seq(2,len,2)]])

  rv <- list(lengths=lengths,
             starts =starts,
             peaks = peaks,
             ends = ends)
  return(rv)

}



##--------------------------------------------------------------------------

area_swe_series <- function(gcm,var.name,swe.file,swe.dir,clip.shp) {

   swe.data <- brick(paste0(swe.dir,swe.file))
   swe.mask <- mask(swe.data,clip.shp)
   swe.series <- as.numeric(cellStats(swe.mask,mean,na.rm=T))*1000
   swe.cover <- as.numeric(cellStats(swe.mask > 0.1,sum,na.rm=T)) ##Greater than 100mm
   swe.volume <- as.numeric(cellStats(swe.mask*800*800,sum,na.rm=T)/10e6) ##In Millions of cubic metres
   nc <- nc_open(paste0(swe.dir,swe.file))
   var.time <- format(netcdf.calendar(nc),'%Y-%m-%d')
   nc_close(nc)

  if (grepl('HadGEM',swe.file)) {
     years <- 1950:2100
     had.years <- as.numeric(format(as.Date(var.time),'%Y'))
     had.ix <- years %in% had.years
     had.fill <- had.cover <- had.volume <- rep(NA,length(years))
     had.fill[had.ix] <- as.numeric(swe.series)
     had.cover[had.ix] <- as.numeric(swe.cover)
     had.volume[had.ix] <- as.numeric(swe.volume)
     had.time <- paste0(years,'-',format(as.Date(var.time[20]),'%m-%d'))
     swe.series <- had.fill
     swe.cover <- had.cover
     swe.volume <- had.volume
     var.time <- had.time
  }
  rv <- list(time=var.time,series=swe.series,cover=swe.cover,volume=swe.volume)
  return(rv)
}



##---------------------------------------------------------------------
read.dir <- '/storage/data/climate/downscale/BCCAQ2+PRISM/bccaq2_tps/BCCAQ2/snow_model/'
gcm.list <- c('ACCESS1-0','CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0','GFDL-ESM2G',
              'HadGEM2-CC','HadGEM2-ES','inmcm4','MIROC5','MPI-ESM-LR','MRI-CGCM3')

shp <- readOGR('/storage/data/projects/rci/data/assessments/metro_van/shapefiles/','MVWaterSheds', stringsAsFactors=F)
clip.shp <- spTransform(shp,CRS("+init=epsg:4326"))
var.name <- 'swe'
type <- 'april_first'

swe.series <- swe.cover <- swe.volume <- matrix(NA,nrow=151,ncol=length(gcm.list))
swe.time <- vector(mode='list',length=length(gcm.list))
temp.years <- read.csv('/storage/data/projects/rci/data/winter_sports/temperature.anomaly.years.csv',header=TRUE,as.is=T)

if (1==1) {

past.cover <- one.cover <- two.cover <- three.cover <- matrix(0,nrow=31,ncol=length(gcm.list))
past.swe <- one.swe <- two.swe <- three.swe <- matrix(0,nrow=31,ncol=length(gcm.list))
past.volume <- one.volume <- two.volume <- three.volume <- matrix(0,nrow=31,ncol=length(gcm.list))

for (i in seq_along(gcm.list)) {
   gcm <- gcm.list[i]
   print(gcm)
   swe.dir <- paste0(read.dir,'calibrated_',gcm,'_PNWNAmet_prism_tps/')
   swe.file <- list.files(path=swe.dir,pattern=type)
   swe.info <-  area_swe_series(gcm,var.name,swe.file,swe.dir,clip.shp)

   swe.series[,i] <- swe.info$series
   swe.cover[,i] <- swe.info$cover
   swe.volume[,i] <- swe.info$volume
   swe.time[[i]] <- swe.info$time   
   swe.years <- format(as.Date(swe.info$time),'%Y')
 
   past.cover[,i] <- swe.info$cover[head(grep(1980,swe.years),1):tail(grep(2010,swe.years),1)]
   past.swe[,i] <- swe.info$series[head(grep(1980,swe.years),1):tail(grep(2010,swe.years),1)]
   past.volume[,i] <- swe.info$volume[head(grep(1980,swe.years),1):tail(grep(2010,swe.years),1)]
   
   ##Intervals for standard temperature changes
   tx <- which(temp.years[,1] == gcm)
      
   one.cover[,i] <- swe.info$cover[head(grep(temp.years[tx,2]-15,swe.years),1):
                              tail(grep(temp.years[tx,2]+15,swe.years),1)]   
   two.cover[,i] <- swe.info$cover[head(grep(temp.years[tx,3]-15,swe.years),1):
                              tail(grep(temp.years[tx,3]+15,swe.years),1)]   
   three.cover[,i] <- swe.info$cover[head(grep(temp.years[tx,4]-15,swe.years),1):
                              tail(grep(temp.years[tx,4]+15,swe.years),1)]   

   one.swe[,i] <- swe.info$series[head(grep(temp.years[tx,2]-15,swe.years),1):
                              tail(grep(temp.years[tx,2]+15,swe.years),1)]   
   two.swe[,i] <- swe.info$series[head(grep(temp.years[tx,3]-15,swe.years),1):
                              tail(grep(temp.years[tx,3]+15,swe.years),1)]   
   three.swe[,i] <- swe.info$series[head(grep(temp.years[tx,4]-15,swe.years),1):
                              tail(grep(temp.years[tx,4]+15,swe.years),1)]     

   one.volume[,i] <- swe.info$volume[head(grep(temp.years[tx,2]-15,swe.years),1):
                              tail(grep(temp.years[tx,2]+15,swe.years),1)]   
   two.volume[,i] <- swe.info$volume[head(grep(temp.years[tx,3]-15,swe.years),1):
                              tail(grep(temp.years[tx,3]+15,swe.years),1)]   
   three.volume[,i] <- swe.info$volume[head(grep(temp.years[tx,4]-15,swe.years),1):
                              tail(grep(temp.years[tx,4]+15,swe.years),1)]     
}

}

##958 cells in the Watershed region


##---------------------------------------------------------------------
##Start, Peak, End, Lengths come from the full time series

reg.dir <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow_series/'

past.starts <- one.starts <- two.starts <- three.starts <- matrix(0,nrow=31,ncol=length(gcm.list))
past.ends <- one.ends <- two.ends <- three.ends <- matrix(0,nrow=31,ncol=length(gcm.list))
past.lengths <- one.lengths <- two.lengths <- three.lengths <- matrix(0,nrow=31,ncol=length(gcm.list))
past.peaks <- one.peaks <- two.peaks <- three.peaks <- matrix(0,nrow=31,ncol=length(gcm.list))

for (i in seq_along(gcm.list)) {
   print(gcm.list[i])
   reg.file <- paste0(reg.dir,'mv_watersheds_',gcm.list[i],'_PNWNAmet_PRISM_TPS_snow_model_data.csv')  
   reg.data <- read.csv(reg.file,header=T,as.is=T)

   ##Scale data 
   years <- format(as.Date(reg.data$Dates),'%Y')
   if (grepl('HadGEM',gcm.list[i])) {
      years[is.na(years)] <- years[which(is.na(years))-2]
   }
   yr.fac <- as.factor(years)
   yr.min <- tapply(reg.data$SWE,yr.fac,min)
   yrs <- unique(years)
   swe.scale <- reg.data$SWE * 0
   for (k in seq_along(yrs)) {
      ix <- grep(yrs[k],years)
      swe.scale[ix] <- reg.data$SWE[ix] - yr.min[k]  
   }

   ##plot(as.Date(reg.data$Dates[10001:23000]),reg.data$SWE[10001:23000],pch=16,main=gcm.list[i],ylim=c(0,2500))
   ##lines(as.Date(reg.data$Dates[10001:23000]),swe.scale[10001:23000],col='red')
   ##abline(h=30)
   
   reg.dates <- get_snow_season_dates(swe.scale,reg.data$Dates)
   starts <- reg.dates$starts
   ends <- reg.dates$ends
   peaks <- reg.dates$peaks
   lengths <- reg.dates$lengths

   past.starts[,i] <- format(starts[starts >= as.Date('1980-01-01') & starts <= as.Date('2010-12-31')],'%j')
   past.ends[,i] <- format(ends[ends >= as.Date('1980-01-01') & ends <= as.Date('2010-12-31')],'%j')
   past.lengths[,i] <- lengths[starts >= as.Date('1980-01-01') & starts <= as.Date('2010-12-31')]
   past.peaks[,i] <- peaks[as.numeric(names(peaks)) >= 1980 & as.numeric(names(peaks)) <= 2010]
   

   ost <- as.Date(paste0(temp.years[i,2]-15,'-01-01'))   
   oen <- as.Date(paste0(temp.years[i,2]+15,'-12-31'))   
   one.starts[,i] <- format(starts[starts > ost & starts <= oen],'%j')
   one.ends[,i] <- format(ends[ends > ost & ends <= oen],'%j')
   one.lengths[,i] <- lengths[starts > ost & starts <= oen]
   one.peaks[,i] <- peaks[as.numeric(names(peaks)) >= temp.years[i,2]-15 &
                          as.numeric(names(peaks)) <= temp.years[i,2]+15]

   tst <- as.Date(paste0(temp.years[i,3]-15,'-01-01'))   
   ten <- as.Date(paste0(temp.years[i,3]+15,'-12-31'))   
   two.starts[,i] <- format(starts[starts > tst & starts <= ten],'%j')
   two.ends[,i] <- format(ends[ends > tst & ends <= ten],'%j')
   two.lengths[,i] <- lengths[starts > tst & starts <= ten]
   two.peaks[,i] <- peaks[as.numeric(names(peaks)) >= temp.years[i,3]-15 &
                          as.numeric(names(peaks)) <= temp.years[i,3]+15]

   est <- as.Date(paste0(temp.years[i,4]-15,'-01-01'))   
   een <- as.Date(paste0(temp.years[i,4]+15,'-12-31'))   
   three.starts[,i] <- c(format(starts[starts > est & starts <= een],'%j'),NA)[1:31]
   three.ends[,i] <- c(format(ends[ends > est & ends <= een],'%j'),NA)[1:31]
   three.lengths[,i] <- c(lengths[starts > est & starts <= een],NA)[1:31]   
   three.peaks[,i] <- c(peaks[as.numeric(names(peaks)) >= temp.years[i,4]-15 &
                            as.numeric(names(peaks)) <= temp.years[i,4]+15],NA)[1:31]
}

past.lengths[past.lengths > 365] <- 365
one.lengths[one.lengths > 365] <- 365
two.lengths[two.lengths > 365] <- 365
three.lengths[three.lengths > 365] <- 365

past.starts <- apply(past.starts,2,as.numeric)
past.ends <- apply(past.ends,2,as.numeric)
one.starts <- apply(one.starts,2,as.numeric)
one.ends <- apply(one.ends,2,as.numeric)
two.starts <- apply(two.starts,2,as.numeric)
two.ends <- apply(two.ends,2,as.numeric)
three.starts <- apply(three.starts,2,as.numeric)
three.ends <- apply(three.ends,2,as.numeric)


##---------------------------------------------------------------------
plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/plots/metro.van.watersheds.boxplots.png')
png(file=plot.file,width=6,height=5,units='in',res=600,pointsize=6,bg='white')

par(mfrow=c(3,2))
par(mar=c(5,6,2,2))
##April 1 Avg SWE

plot(c(),xlim=c(0.5,4.5),ylim=c(0,1700),axes=FALSE,
         xlab='Temperature Change',ylab='April 1 SWE (mm)',
         cex.lab=1.95,cex.axis=1.95,yaxs='i')
abline(h=seq(250,2000,250),col='lightgray',lty=2)
boxplot(apply(past.swe,2,mean,na.rm=T),at=1,col='blue',add=T,axes=F)
boxplot(apply(one.swe,2,mean,na.rm=T),at=2,add=T,col='goldenrod',axes=F)
boxplot(apply(two.swe,2,mean,na.rm=T),at=3,add=T,col='orange',axes=F)
boxplot(apply(three.swe,2,mean,na.rm=T),at=4,add=T,col='red',axes=F)

axis(1,at=c(1,2,3,4),label=c('Historic','1\u00B0C','2\u00B0C','3\u00B0C'),cex.axis=1.95,cex=1.95,cex.lab=1.95)
axis(2,at=seq(0,1750,250),label=seq(0,1750,250),cex.axis=1.95)
box(which='plot',lwd=1.5)

##April 1 Snow Cover
plot(c(),xlim=c(0.5,4.5),ylim=c(30,90),axes=FALSE,
         xlab='Temperature Change',ylab='Arpil 1 Snow Cover (%)',
         cex.lab=1.95,cex.axis=1.95,yaxs='i')
abline(h=seq(40,100,20),col='lightgray',lty=2)
boxplot(apply(past.cover,2,mean,na.rm=T)/958*100,at=1,col='blue',add=T,axes=F)
boxplot(apply(one.cover,2,mean,na.rm=T)/958*100,at=2,add=T,col='goldenrod',axes=F)
boxplot(apply(two.cover,2,mean,na.rm=T)/958*100,at=3,add=T,col='orange',axes=F)
boxplot(apply(three.cover,2,mean,na.rm=T)/958*100,at=4,add=T,col='red',axes=F)
axis(1,at=c(1,2,3,4),label=c('Historic','1\u00B0C','2\u00B0C','3\u00B0C'),cex.axis=1.95,cex=1.95,cex.lab=1.95)
axis(2,at=seq(0,100,20),label=seq(0,100,20),cex.axis=1.95)
box(which='plot',lwd=1.5)

##Starts
plot(c(),xlim=c(0.5,4.5),ylim=c(274,325),axes=FALSE,
         xlab='Temperature Change',ylab='Start Day',
         cex.lab=1.95,cex.axis=1.95,yaxs='i')
abline(h=c(274,288,305,319),col='lightgray',lty=2)
boxplot(apply(past.starts,2,mean,na.rm=T),at=1,col='blue',add=T,axes=F)
boxplot(apply(one.starts,2,mean,na.rm=T),at=2,add=T,col='goldenrod',axes=F)
boxplot(apply(two.starts,2,mean,na.rm=T),at=3,add=T,col='orange',axes=F)
boxplot(apply(three.starts,2,mean,na.rm=T),at=4,add=T,col='red',axes=F)
axis(1,at=c(1,2,3,4),label=c('Historic','1\u00B0C','2\u00B0C','3\u00B0C'),cex.axis=1.95,cex=1.95,cex.lab=1.95)
##axis(2,at=seq(275,325,15),label=seq(275,325,15),cex.axis=1.95)
axis(2,at=c(274,288,305,319),label=c('Oct 1','Oct 15','Nov 1','Nov 15'),cex.axis=1.95)
box(which='plot',lwd=1.5)

##Peaks
plot(c(),xlim=c(0.5,4.5),ylim=c(32,115),axes=FALSE,
         xlab='Temperature Change',ylab='Peak SWE Day',
         cex.lab=1.95,cex.axis=1.95,yaxs='i')
abline(h=c(32,60,91,121),col='lightgray',lty=2)
boxplot(apply(past.peaks,2,mean,na.rm=T),at=1,col='blue',add=T,axes=F)
boxplot(apply(one.peaks,2,mean,na.rm=T),at=2,add=T,col='goldenrod',axes=F)
boxplot(apply(two.peaks,2,mean,na.rm=T),at=3,add=T,col='orange',axes=F)
boxplot(apply(three.peaks,2,mean,na.rm=T),at=4,add=T,col='red',axes=F)
axis(1,at=c(1,2,3,4),label=c('Historic','1\u00B0C','2\u00B0C','3\u00B0C'),cex.axis=1.95,cex=1.95,cex.lab=1.95)
axis(2,at=c(32,60,91),label=c('Feb','Mar','Apr'),cex.axis=1.95)
##axis(2,at=seq(25,125,15),label=seq(25,125,15),cex.axis=1.95)
box(which='plot',lwd=1.5)

##Ends
plot(c(),xlim=c(0.5,4.5),ylim=c(182,290),axes=FALSE,
         xlab='Temperature Change',ylab='End Date',
         cex.lab=1.95,cex.axis=1.95,yaxs='i')
abline(h=c(152,182,213,244,274),col='lightgray',lty=2)
boxplot(apply(past.ends,2,mean,na.rm=T),at=1,col='blue',add=T,axes=F)
boxplot(apply(one.ends,2,mean,na.rm=T),at=2,add=T,col='goldenrod',axes=F)
boxplot(apply(two.ends,2,mean,na.rm=T),at=3,add=T,col='orange',axes=F)
boxplot(apply(three.ends,2,mean,na.rm=T),at=4,add=T,col='red',axes=F)
axis(1,at=c(1,2,3,4),label=c('Historic','1\u00B0C','2\u00B0C','3\u00B0C'),cex.axis=1.95,cex=1.95,cex.lab=1.95)
##axis(2,at=seq(150,300,25),label=seq(150,300,25),cex.axis=1.95)
axis(2,at=c(152,182,213,244,274) ,label=c('Jun','Jul','Aug','Sep','Oct'),cex.axis=1.95)
box(which='plot',lwd=1.5)

##Lengths
plot(c(),xlim=c(0.5,4.5),ylim=c(250,370),axes=FALSE,
         xlab='Temperature Change',ylab='Length (Days)',
         cex.lab=1.95,cex.axis=1.95,yaxs='i')
abline(h=seq(225,350,25),col='lightgray',lty=2)
boxplot(apply(past.lengths,2,mean,na.rm=T),at=1,col='blue',add=T,axes=F)
boxplot(apply(one.lengths,2,mean,na.rm=T),at=2,add=T,col='goldenrod',axes=F)
boxplot(apply(two.lengths,2,mean,na.rm=T),at=3,add=T,col='orange',axes=F)
boxplot(apply(three.lengths,2,mean,na.rm=T),at=4,add=T,col='red',axes=F)
axis(1,at=c(1,2,3,4),label=c('Historic','1\u00B0C','2\u00B0C','3\u00B0C'),cex.axis=1.95,cex=1.95,cex.lab=1.95)
axis(2,at=seq(225,350,25),label=seq(225,350,25),cex.axis=1.95)
box(which='plot',lwd=1.5)


dev.off()
