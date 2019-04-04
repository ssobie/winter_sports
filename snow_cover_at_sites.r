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
snodas.snow.days <- vector(mode='list',length=slen)

era.snodas.days <- vector(mode='list',length=slen)
ncep2.snodas.days <- vector(mode='list',length=slen)
modis.snodas.days <- vector(mode='list',length=slen)


elevs <- rep(0,slen)
lons <- rep(0,slen)
##---------------------------------------------------------------------------------
##MODIS 
modis.dir <- '/storage/data/projects/rci/data/winter_sports/MODIS_VAN_WHISTLER/'
mnc.file <- paste0(modis.dir,'snc.modis.van_whistler.20010101-20181231.nc')
mnc.nc <- nc_open(mnc.file)
m.lon <- ncvar_get(mnc.nc,'lon')
m.lat <- ncvar_get(mnc.nc,'lat')
modis.time <- netcdf.calendar(mnc.nc)

##---------------------------------------------------------------------------------
##SNODAS 
snodas.dir <- '/storage/data/projects/rci/data/winter_sports/obs/SNODAS/ncdf4_files/'
snc.file <- paste0(snodas.dir,'swe_snodas_modis_grid_van_whistler_20100101-20181231.nc')
snc.nc <- nc_open(snc.file)
s.lon <- ncvar_get(snc.nc,'lon')
s.lat <- ncvar_get(snc.nc,'lat')
snodas.time <- netcdf.calendar(snc.nc)

##---------------------------------------------------------------------------------
##SNOW MODEL from Simulations
snow.time.dir <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/')
esnw.time.file <- paste0(snow.time.dir,'snowdepth_BCCAQ2-PRISM_ERA_19790101-20181031.nc')
esnw.nc <- nc_open(esnw.time.file)
era.time <- netcdf.calendar(esnw.nc)

nsnw.time.file <- paste0(snow.time.dir,'snowdepth_BCCAQ2-PRISM_NCEP2_19790101-20181031.nc')
nsnw.nc <- nc_open(nsnw.time.file)
ncep2.time <- netcdf.calendar(nsnw.nc)

ncep2.match <- format(ncep2.time,'%Y-%m-%d') %in% format(era.time,'%Y-%m-%d')


##---------------------------------------------------------------------------------

date.match <- format(era.time,'%Y-%m-%d') %in% format(modis.time,'%Y-%m-%d')
modis.match <- format(modis.time,'%Y-%m-%d') %in% format(era.time,'%Y-%m-%d')
snodas.match <- format(modis.time,'%Y-%m-%d') %in% format(snodas.time,'%Y-%m-%d')

era.snodas.match <- format(era.time,'%Y-%m-%d') %in% format(snodas.time,'%Y-%m-%d')
snodas.era.match <- format(snodas.time,'%Y-%m-%d') %in% format(era.time,'%Y-%m-%d')

modis.era.snodas.match <- format(modis.time,'%Y-%m-%d') %in% format(era.time[era.snodas.match],'%Y-%m-%d')


bnds <- range(which(date.match))

for (i in 1:slen) {
  site <- sites[i]
  snow.dir <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sims/')
  era.file <- paste0(snow.dir,site,'.era.snow.1001.csv')
  era.sims <- as.matrix(read.csv(era.file,header=TRUE,as.is=TRUE))

  ncep2.file <- paste0(snow.dir,site,'.ncep2.snow.1001.csv')
  ncep2.sims <- as.matrix(read.csv(ncep2.file,header=TRUE,as.is=TRUE))[ncep2.match,]


  print(site)
  coords <- get.coordinates(site)
  lons[i] <- coords[1]
  elevs[i] <- coords[3]
  
  m.ox <- which.min(abs(coords[1]-m.lon))
  m.ax <- which.min(abs(coords[2]-m.lat))

  mnc.data <- ncvar_get(mnc.nc,'snc',start=c(m.ox,m.ax,1),count=c(1,1,-1))[modis.match]
  mnc.modis <- mnc.data
  mnc.modis[mnc.modis > 0 & mnc.modis <=100] <- 1
  mnc.modis[mnc.modis > 100] <- NA
  mnc.valid <- mnc.modis
  ##mnc.modis[mnc.modis == 0]  <- NA

  s.ox <- which.min(abs(coords[1]-s.lon))
  s.ax <- which.min(abs(coords[2]-s.lat))

  snc.data <- ncvar_get(snc.nc,'swe',start=c(s.ox,s.ax,1),count=c(1,1,-1))
  snc.snodas <- snc.data
  snc.snodas[snc.snodas <= 5] <- 0
  snc.snodas[snc.snodas > 5] <- 1

  era.cover <- era.sims ##[bnds[1]:bnds[2],]
  era.cover[era.cover>0.05] <- 1
  era.cover[era.cover<=0.05] <- 0

  ncep2.cover <- ncep2.sims ##[bnds[1]:bnds[2],]
  ncep2.cover[ncep2.cover>0.05] <- 1
  ncep2.cover[ncep2.cover<=0.05] <- 0
 
  mnc.matrix <- matrix(mnc.modis,nrow=length(mnc.modis),ncol=1001)         
  snc.matrix <- matrix(snc.snodas[snodas.era.match],nrow=sum(snodas.era.match),ncol=1001)         
    
  ##era.diff <- era.cover - mnc.matrix
  ##ncep2.diff <- ncep2.cover - mnc.matrix
  ##era.cover.diff[[i]] <- apply(era.diff,2,function(x){round(sum(x==0,na.rm=T)/sum(!is.na(x))*100,1)})
  ##ncep2.cover.diff[[i]] <- apply(ncep2.diff,2,function(x){round(sum(x==0,na.rm=T)/sum(!is.na(x))*100,1)})

  ##For snow only

  flag <- is.na(mnc.matrix)
  era.snow <- era.cover[bnds[1]:bnds[2],]
  ncep2.snow <- ncep2.cover[bnds[1]:bnds[2],] 
  era.snow[flag] <- NA
  ncep2.snow[flag] <- NA

  mflag <- is.na(mnc.modis[snodas.match])
  modis.snodas <- mnc.modis[snodas.match]
  modis.era.snodas <- mnc.modis[modis.era.snodas.match]
  sflag <- is.na(modis.era.snodas)

  era.snodas <- era.cover[era.snodas.match,]
  era.snodas[sflag,] <- NA
  ncep2.snodas <- ncep2.cover[era.snodas.match,]
  ncep2.snodas[sflag,] <- NA  

  era.cover.diff[[i]] <- (apply(era.snow,2,sum,na.rm=T) - sum(mnc.modis,na.rm=T))/sum(mnc.modis,na.rm=T)*100
  ncep2.cover.diff[[i]] <- (apply(ncep2.snow,2,sum,na.rm=T) - sum(mnc.modis,na.rm=T))/sum(mnc.modis,na.rm=T)*100                            

  era.snow.days[[i]] <- apply(era.snow,2,sum,na.rm=T)
  ncep2.snow.days[[i]] <- apply(ncep2.snow,2,sum,na.rm=T)  
  modis.snow.days[[i]] <- sum(mnc.modis,na.rm=T)
  modis.valid.days[[i]] <- sum(!is.na(mnc.modis))


  era.snodas.days[[i]] <- apply(era.snodas,2,sum,na.rm=T)
  ncep2.snodas.days[[i]] <- apply(ncep2.snodas,2,sum,na.rm=T)  
  modis.snodas.days[[i]] <- sum(modis.snodas,na.rm=T)  
  snodas.snow.days[[i]] <- sum(snc.snodas[!mflag],na.rm=T)


  print(paste0('Total MODIS Observable Days: ',sum(!is.na(mnc.modis))))
  print(paste0('Total MODIS Snow Days: ',sum(mnc.modis,na.rm=T)))
  print(paste0('Total SNODAS Snow Days: ',sum(snc.snodas,na.rm=T)))
  print(paste0('Average ERA Snow Days: ',mean(apply(era.snow,2,sum,na.rm=T))))
  print(paste0('Average NCEP2 Snow Days: ',mean(apply(ncep2.snow,2,sum,na.rm=T))))


}

nc_close(mnc.nc)
nc_close(snc.nc)
nc_close(nsnw.nc)
nc_close(esnw.nc)


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
if (1==1) {
plot.file <- paste0('/storage/data/projects/rci/data/winter_sports/plots/',model,'.total.snow.days.by.longitude.sites.1001.2018.png')
plot.title <- '' ##paste0('Total Snow Days')
leg.title <- 'Days'
ranked.elevs <- order(elevs)
ranked.lons <- order(lons)
png(filename=plot.file,width=1000,height=700)
layout(mat = matrix(c(1,2),
                    nrow = 2,
                    ncol = 1),
             heights = c(1, 1.5),    # Heights of the two rows
             widths = 1)     # Widths of the two columns

par(mar=c(0,5,3,3))
plot(0:slen,0:slen,xlab='',ylab='Snow Days',yaxs='i',
     col='white',main=plot.title,cex.axis=1.75,cex.lab=1.75,cex.main=2,
     xlim=c(1,slen),ylim=c(0,800),axes=FALSE)        
##axis(1,at=1:slen,site.names[ranked.lons],cex=1.75,cex.axis=1.75,las=2)                                                        
axis(2,at=seq(0,1200,100),seq(0,1200,100),cex=1.75,cex.axis=1.75)                                                        
abline(h=seq(0,1200,100),lty=2,col='gray',lwd=2)
abline(v=1:slen,col='gray')
for (j in 1:slen) {
    print(elevs[ranked.lons[j]])
    boxplot(at=j-0.175,x=era.snodas.days[[ranked.lons[j]]],add=TRUE,axes=F,boxwex=0.7,col='blue',border='blue')
    boxplot(at=j+0.175,x=ncep2.snodas.days[[ranked.lons[j]]],add=TRUE,axes=F,boxwex=0.7,col='green',border='green')
    points(x=j,y=modis.snodas.days[[ranked.lons[j]]],pch='-',cex=5,col='black')
    points(x=j,y=snodas.snow.days[[ranked.lons[j]]],pch='-',cex=5,col='red')
}
rect(0,710,2.3,800,border='black',col='white')
text(x=1.3,y=750,'2009-2018',cex=1.5)

##legend('topleft',leg=c('MODIS','SNODAS','ERA-I','NCEP2'),col=c('black','red','blue','green'),pch=15,cex=1.5)
box(which='plot')

par(mar=c(10,5,0.1,3))
plot(0:slen,0:slen,xlab='',ylab='Snow Days',yaxs='i',
     col='white',main=plot.title,cex.axis=1.75,cex.lab=1.75,cex.main=2,
     xlim=c(1,slen),ylim=c(200,1100),axes=FALSE)        
axis(1,at=1:slen,site.names[ranked.lons],cex=1.75,cex.axis=1.75,las=2)                                                        
axis(2,at=c(300,600,900),c(300,600,900),cex=1.75,cex.axis=1.75)                                                        
abline(h=seq(0,1200,100),lty=2,col='gray',lwd=2)
abline(v=1:slen,col='gray')
for (j in 1:slen) {
    print(elevs[ranked.lons[j]])
    boxplot(at=j-0.175,x=era.snow.days[[ranked.lons[j]]],add=TRUE,axes=F,boxwex=0.7,col='blue',border='blue')
    boxplot(at=j+0.175,x=ncep2.snow.days[[ranked.lons[j]]],add=TRUE,axes=F,boxwex=0.7,col='green',border='green')
    points(x=j,y=modis.snow.days[[ranked.lons[j]]],pch='-',cex=5,col='black')
}
rect(0,1010,2.3,1100,border='black',col='white')
text(x=1.3,y=1050,'2001-2018',cex=1.5)
##legend('topleft',leg=c('MODIS','ERA-I','NCEP2'),col=c('black','blue','green'),pch=15,cex=1.5)
legend('bottomleft',leg=c('MODIS','SNODAS','ERA-I','NCEP2'),col=c('black','red','blue','green'),pch=15,cex=1.5)
box(which='plot')
dev.off()

}

##MODIS Total Snow Days
if (1==0) {
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
