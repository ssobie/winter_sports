##Script to plot the fraction of useable MODIS data
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)
source('/home/ssobie/code/repos/winter_sports/site.coordinates.r',chdir=T)

library(ncdf4)
library(PCICt)

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
           'upper_squamish',
           'spuzzum_creek',
           'chilliwack_river',
           'tenquille_lake')

site.names <- c('Shovelnose\nMountain',
                'Brookmere', 
                'Lightning\nLake',
                'Callaghan',
                'Orchid\nLake',
                'Palisade\nLake',
                'Grouse\nMountain',
                'Dog\nMountain',
                'Stave\nLake',
                'Nahatlatch',
                'Wahleach',
                'Klesilkwa',
                'Hamilton\nHill',
                'Upper\nSquamish',
                'Spuzzum\nCreek',
                'Chilliwack\nRiver',
                'Tenquille\nLake')


##sites <- 'klesilkwa'
slen <- length(sites)

model <- 'ERA'

era.cover.diff <- vector(mode='list',length=slen)
ncep2.cover.diff <- vector(mode='list',length=slen)

era.snow.days <- vector(mode='list',length=slen)
ncep2.snow.days <- vector(mode='list',length=slen)
modis.snow.days <- vector(mode='list',length=slen)
modis.valid.days <- vector(mode='list',length=slen)

elevs <- rep(0,slen)
##---------------------------------------------------------------------------------
##MODIS 
modis.dir <- '/storage/data/projects/rci/data/winter_sports/MODIS_VAN_WHISTLER/'
snc.file <- paste0(modis.dir,'snc.modis.van_whistler.20010101-20151231.nc')
snc.nc <- nc_open(snc.file)
m.lon <- ncvar_get(snc.nc,'lon')
m.lat <- ncvar_get(snc.nc,'lat')
modis.time <- netcdf.calendar(snc.nc)

##---------------------------------------------------------------------------------
##SNOW MODEL from Simulations
snow.time.dir <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/',model,'/')
snw.time.file <- paste0(snow.time.dir,'snowdepth_BCCAQ-PRISM_',model,'_rcp85_r1_1979-2016.nc')
snw.nc <- nc_open(snw.time.file)
s.lon <- ncvar_get(snw.nc,'lon')
s.lat <- ncvar_get(snw.nc,'lat')
snow.time <- netcdf.calendar(snw.nc)

##---------------------------------------------------------------------------------

date.match <- format(snow.time,'%Y-%m-%d') %in% format(modis.time,'%Y-%m-%d')
bnds <- range(which(date.match))

for (i in 1:slen) {
  site <- sites[i]
  snow.dir <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/snow_sims/')
  era.file <- paste0(snow.dir,site,'.era.snow.11.csv')
  era.sims <- read.csv(era.file,header=TRUE,as.is=TRUE)

  ncep2.file <- paste0(snow.dir,site,'.ncep2.snow.11.csv')
  ncep2.sims <- read.csv(ncep2.file,header=TRUE,as.is=TRUE)

  print(site)
  coords <- get.coordinates(site)
  elevs[i] <- coords[3]
  m.ox <- which.min(abs(coords[1]-m.lon))
  m.ax <- which.min(abs(coords[2]-m.lat))

  s.ox <- which.min(abs(coords[1]-s.lon))
  s.ax <- which.min(abs(coords[2]-s.lat))

  snc.data <- ncvar_get(snc.nc,'snc',start=c(m.ox,m.ax,1),count=c(1,1,-1))
  snc.modis <- snc.data
  snc.modis[snc.modis > 0 & snc.modis <=100] <- 1
  snc.modis[snc.modis > 100] <- NA
  snc.valid <- snc.modis
  ##snc.modis[snc.modis == 0]  <- NA

  ##snw.data <- ncvar_get(snw.nc,'snowdepth',start=c(s.ox,s.ax,bnds[1]),count=c(1,1,diff(bnds)+1))
  era.cover <- era.sims[bnds[1]:bnds[2],]
  era.cover[era.cover>0.05] <- 1
  era.cover[era.cover<=0.05] <- 0

  ncep2.cover <- ncep2.sims[bnds[1]:bnds[2],]
  ncep2.cover[ncep2.cover>0.05] <- 1
  ncep2.cover[ncep2.cover<=0.05] <- 0
 
  snc.matrix <- matrix(snc.modis,nrow=length(snc.modis),ncol=11)         
  ##era.diff <- era.cover - snc.matrix
  ##ncep2.diff <- ncep2.cover - snc.matrix
  ##era.cover.diff[[i]] <- apply(era.diff,2,function(x){round(sum(x==0,na.rm=T)/sum(!is.na(x))*100,1)})
  ##ncep2.cover.diff[[i]] <- apply(ncep2.diff,2,function(x){round(sum(x==0,na.rm=T)/sum(!is.na(x))*100,1)})

  ##For snow only
  flag <- is.na(snc.matrix)
  era.snow <- era.cover
  ncep2.snow <- ncep2.cover 
  era.snow[flag] <- NA
  ncep2.snow[flag] <- NA

  era.cover.diff[[i]] <- (apply(era.snow,2,sum,na.rm=T) - sum(snc.modis,na.rm=T))/sum(snc.modis,na.rm=T)*100
  ncep2.cover.diff[[i]] <- (apply(ncep2.snow,2,sum,na.rm=T) - sum(snc.modis,na.rm=T))/sum(snc.modis,na.rm=T)*100                            

  era.snow.days[[i]] <- apply(era.snow,2,sum,na.rm=T)
  ncep2.snow.days[[i]] <- apply(ncep2.snow,2,sum,na.rm=T)  
  modis.snow.days[[i]] <- sum(snc.modis,na.rm=T)
  modis.valid.days[[i]] <- sum(!is.na(snc.modis))

  print(paste0('Total MODIS Observable Days: ',sum(!is.na(snc.modis))))
  print(paste0('Total MODIS Snow Days: ',sum(snc.modis,na.rm=T)))
  print(paste0('Average ERA Snow Days: ',mean(apply(era.snow,2,sum,na.rm=T))))
  print(paste0('Average NCEP2 Snow Days: ',mean(apply(ncep2.snow,2,sum,na.rm=T))))


}

nc_close(snc.nc)
nc_close(snw.nc)


##Model Total Snow Cover
if (1==0) {
plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/plots/',model,'.snow.only.prct.at.sites.11.png')
plot.title <- paste0('MODIS Snow Only Comparison')
leg.title <- 'Percent (%)'
ranked.elevs <- order(elevs)
png(filename=plot.file,width=1000,height=500)
par(mar=c(10,5,3,3))
plot(0:slen,0:slen,xlab='',ylab='Snow Cover Match (%)',yaxs='i',
     col='white',main=plot.title,cex.axis=1.75,cex.lab=1.75,cex.main=2,
     xlim=c(1,slen),ylim=c(-20,70),axes=FALSE)        
axis(1,at=1:slen,site.names[ranked.elevs],cex=1.75,cex.axis=1.75,las=2)                                                        
axis(2,at=seq(-20,70,20),seq(-20,70,20),cex=1.75,cex.axis=1.75)                                                        
abline(h=seq(-20,70,20),lty=2,col='gray',lwd=2)
abline(v=1:slen,col='gray')
for (j in 1:slen) {
    print(elevs[ranked.elevs[j]])
    boxplot(at=j-0.15,x=era.cover.diff[[ranked.elevs[j]]],add=TRUE,axes=F,boxwex=0.6)
    boxplot(at=j+0.15,x=ncep2.cover.diff[[ranked.elevs[j]]],add=TRUE,axes=F,boxwex=0.6)
}

box(which='plot')
dev.off()

}

##MODIS Total Snow Days
if (1==0) {
plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/plots/',model,'.total.snow.days.prct.at.sites.11.png')
plot.title <- paste0('Total Snow Days')
leg.title <- 'Days'
ranked.elevs <- order(elevs)
png(filename=plot.file,width=1000,height=500)
par(mar=c(10,5,3,3))
plot(0:slen,0:slen,xlab='',ylab='Snow Days',yaxs='i',
     col='white',main=plot.title,cex.axis=1.75,cex.lab=1.75,cex.main=2,
     xlim=c(1,slen),ylim=c(200,900),axes=FALSE)        
axis(1,at=1:slen,site.names[ranked.elevs],cex=1.75,cex.axis=1.75,las=2)                                                        
axis(2,at=seq(200,900,100),seq(200,900,100),cex=1.75,cex.axis=1.75)                                                        
abline(h=seq(200,900,100),lty=2,col='gray',lwd=2)
abline(v=1:slen,col='gray')
for (j in 1:slen) {
    print(elevs[ranked.elevs[j]])
    boxplot(at=j-0.15,x=era.snow.days[[ranked.elevs[j]]],add=TRUE,axes=F,boxwex=0.6)
    boxplot(at=j+0.15,x=ncep2.snow.days[[ranked.elevs[j]]],add=TRUE,axes=F,boxwex=0.6)
    points(x=j,y=modis.snow.days[[ranked.elevs[j]]],pch='-',cex=2,col='red')
}

box(which='plot')
dev.off()

}

##MODIS Total Snow Days
if (1==1) {
plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/plots/',model,'.valid.days.at.sites.11.png')
plot.title <- paste0('MODIS Observable Days')
leg.title <- 'Days'
ranked.elevs <- order(elevs)
png(filename=plot.file,width=1000,height=500)
par(mar=c(10,5,3,3))
plot(0:slen,0:slen,xlab='',ylab='Snow Days',yaxs='i',
     col='white',main=plot.title,cex.axis=1.75,cex.lab=1.75,cex.main=2,
     xlim=c(1,slen),ylim=c(1200,2200),axes=FALSE)        
axis(1,at=1:slen,site.names[ranked.elevs],cex=1.75,cex.axis=1.75,las=2)                                                        
axis(2,at=seq(1200,2200,200),seq(1200,2200,200),cex=1.75,cex.axis=1.75)                                                        
abline(h=seq(1200,2200,200),lty=2,col='gray',lwd=2)
abline(v=1:slen,col='gray')
for (j in 1:slen) {
    print(elevs[ranked.elevs[j]])
    points(x=j,y=modis.valid.days[[ranked.elevs[j]]],pch='-',cex=3)
}

box(which='plot')
dev.off()

}
