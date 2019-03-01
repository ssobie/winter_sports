##Script to compare snowdepth, SWE and snow cover between modelled values, 
##snow course sites and MODIS snow cover
library(ncdf4)
library(plotrix)

##source('/storage/home/ssobie/code/repos/winter_sports/test.snow.model.r',chdir=T)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

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
           'hamilton_hill')
ilen <- length(sites)
asps <- c('upper_squamish',
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



course.elevs <- rep(0,length(sites))
pillow.elevs <- rep(0,length(asps))

course.site.swe <- vector(mode='list',length=length(sites))
pillow.site.swe <- vector(mode='list',length=length(asps))
model.site.swe <- vector(mode='list',length=length(sites))

slen <- 1001

era.file <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/grouse_mountain_ERA_800m_data.csv'
era.data <- read.csv(era.file,header=T,as.is=T)
era.dates <- era.data$Dates
era.swe.sims <- matrix(0,nrow=slen,ncol=dim(era.data)[1])

ncep2.file <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/grouse_mountain_NCEP2_800m_data.csv'
ncep2.data <- read.csv(ncep2.file,header=T,as.is=T)
ncep2.dates <- ncep2.data$Dates
ncep2.swe.sims <- matrix(0,nrow=slen,ncol=dim(era.data)[1])

model.match <- format(as.Date(ncep2.dates),'%Y-%m-%d') %in% format(as.Date(era.dates),'%Y-%m-%d') 

ncep2.course.mae <- era.course.mae <- matrix(0,nrow=length(sites),ncol=slen)
ncep2.pillow.mae <- era.pillow.mae <- matrix(0,nrow=length(asps),ncol=slen)
ncep2.snodas.mae <- era.snodas.mae <- matrix(0,nrow=length(sites)+length(asps),ncol=slen)
snodas.site.mae <- rep(0,length(sites))
snodas.asp.mae <- rep(0,length(asps))

model.dir <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sims/'

snodas.file <- "/storage/data/projects/rci/data/winter_sports/obs/SNODAS/ncdf4_files/swe_snodas_modis_grid_van_whistler_20100101-20181231.nc"
snc <- nc_open(snodas.file)
lon <- ncvar_get(snc,'lon')
lat <- ncvar_get(snc,'lat')
snodas.dates <- as.character(netcdf.calendar(snc))

ncep2.snodas.match <-  format(as.Date(ncep2.dates),'%Y-%m-%d') %in% format(as.Date(snodas.dates),'%Y-%m-%d') 
snodas.ncep2.match <-  format(as.Date(snodas.dates),'%Y-%m-%d') %in% format(as.Date(ncep2.dates),'%Y-%m-%d') 

era.snodas.match <-  format(as.Date(era.dates),'%Y-%m-%d') %in% format(as.Date(snodas.dates),'%Y-%m-%d') 
snodas.era.match <-  format(as.Date(snodas.dates),'%Y-%m-%d') %in% format(as.Date(era.dates),'%Y-%m-%d') 


for (i in seq_along(sites)) {
    site <- sites[i]
 
    print(site)
    ##Reanalysis 800m data
    era.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/',site,'_ERA_800m_data.csv')
    era.data <- read.csv(era.file,header=T,as.is=T)
    dates <- era.data$Dates

    coords <- get.coordinates(site)
    lat.bnds <- coords[2]
    course.elevs[i] <- coords[3]

    ##Snow Course Data
    course.file <- paste('/storage/data/projects/rci/data/winter_sports/obs/snow_courses/',site,'_snow_course.csv',sep='')
    course.data <- read.csv(course.file,header=T,as.is=T)
    course.dates <- format(as.Date(course.data[,1]),'%Y-%m-%d')
    course.swe <- course.data[,3] ##mm
    date.subset <- format(as.Date(dates),'%Y-%m-%d') %in% course.dates
    course.subset <- course.dates %in% format(as.Date(dates),'%Y-%m-%d')

    ##SNODAS Data at Courses
    lon.ix <- which.min(abs(coords[1]-lon))
    lat.ix <- which.min(abs(coords[2]-lat))
    snodas.swe <- ncvar_get(snc,'swe',start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
    snodas.date.subset <- format(as.Date(snodas.dates),'%Y-%m-%d') %in% course.dates
    snodas.course.subset <- course.dates %in% format(as.Date(snodas.dates),'%Y-%m-%d')
    snodas.site.mae[i] <- mean(abs(course.swe[snodas.course.subset] - snodas.swe[snodas.date.subset]),na.rm=T)
      
    swe.era <- read.csv(paste0(model.dir,site,'.era.swe.1001.csv'),header=T,as.is=T)
    era.swe.mean <- apply(swe.era,1,mean,na.rm=T)
    swe.ncep2 <- read.csv(paste0(model.dir,site,'.ncep2.swe.1001.csv'),header=T,as.is=T)[model.match,]
    ncep2.swe.mean <- apply(swe.ncep2,1,mean,na.rm=T)

    print(length(course.swe[course.subset]))
    print(dim(swe.era[date.subset,]))
    print(dim(swe.ncep2[date.subset,]))

    for (j in 1:slen) {
      era.course.mae[i,j] <- mean(abs(course.swe[course.subset]-swe.era[date.subset,j]),na.rm=T)
      ncep2.course.mae[i,j] <- mean(abs(course.swe[course.subset]-swe.ncep2[date.subset,j]),na.rm=T)
      era.snodas.mae[i,j] <- mean(abs(snodas.swe[snodas.era.match]-swe.era[era.snodas.match,j]),na.rm=T)      
      ncep2.snodas.mae[i,j] <- mean(abs(snodas.swe[snodas.ncep2.match]-swe.ncep2[ncep2.snodas.match,j]),na.rm=T)      

    }  
}        

for (i in seq_along(asps)) {
    asp <- asps[i]
 
    ##Reanalysis 800m data
    clim.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/',asp,'_ERA_800m_data.csv')
    clim.data <- read.csv(clim.file,header=T,as.is=T)
    dates <- clim.data$Dates

    coords <- get.coordinates(asp)
    lat.bnds <- coords[2]
    pillow.elevs[i] <- coords[3]

    ##Snow Pillow Data
    pillow.file <- paste('/storage/data/projects/rci/data/winter_sports/obs/snow_pillow/',asp,'_asp.csv',sep='')
    pillow.data <- read.csv(pillow.file,header=T,as.is=T)
    pillow.dates <- format(as.Date(pillow.data[,2]),'%Y-%m-%d')
    pillow.swe <- pillow.data[,11] ##mm
    date.subset <- format(as.Date(dates),'%Y-%m-%d') %in% pillow.dates
    pillow.subset <- pillow.dates %in% format(as.Date(dates),'%Y-%m-%d')

    ##SNODAS Data at Courses
    lon.ix <- which.min(abs(coords[1]-lon))
    lat.ix <- which.min(abs(coords[2]-lat))
    snodas.swe <- ncvar_get(snc,'swe',start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
    snodas.date.subset <- format(as.Date(snodas.dates),'%Y-%m-%d') %in% pillow.dates
    snodas.pillow.subset <- pillow.dates %in% format(as.Date(snodas.dates),'%Y-%m-%d')

    swe.era <- read.csv(paste0(model.dir,asp,'.era.swe.1001.csv'),header=T,as.is=T)
    era.swe.mean <- apply(swe.era,1,mean,na.rm=T)
    swe.ncep2 <- read.csv(paste0(model.dir,asp,'.ncep2.swe.1001.csv'),header=T,as.is=T)
    ncep2.swe.mean <- apply(swe.ncep2,1,mean,na.rm=T)

    print(length(pillow.swe[pillow.subset]))
    print(dim(swe.era[date.subset,]))
    print(dim(swe.ncep2[date.subset,]))
 
    snodas.asp.mae[i] <- mean(abs(pillow.swe[snodas.pillow.subset] - snodas.swe[snodas.date.subset]),na.rm=T)      
    for (j in 1:slen) {
      era.pillow.mae[i,j] <- mean(abs(pillow.swe[pillow.subset]-swe.era[date.subset,j]),na.rm=T)
      ncep2.pillow.mae[i,j] <- mean(abs(pillow.swe[pillow.subset]-swe.ncep2[date.subset,j]),na.rm=T)
      era.snodas.mae[i+ilen,j] <- mean(abs(snodas.swe[snodas.era.match]-swe.era[era.snodas.match,j]),na.rm=T)      
      ncep2.snodas.mae[i+ilen,j] <- mean(abs(snodas.swe[snodas.ncep2.match]-swe.ncep2[ncep2.snodas.match,j]),na.rm=T)      
    }  
}        

nc_close(snc)

all.era.mae <- rbind(era.course.mae,era.pillow.mae)
all.ncep2.mae <- rbind(ncep2.course.mae,ncep2.pillow.mae)
all.snodas.mae <- c(snodas.site.mae,snodas.asp.mae)

elevs <- c(course.elevs,pillow.elevs)
alen <- length(sites)+length(asps)

##-------------------------------------------------------------------


##Loop over sites
plot.dir <- '/storage/data/projects/rci/data/winter_sports/plots/'

leg.title <- 'mm'
ranked.elevs <- order(elevs)
png(filename=paste0(plot.dir,'era.ncep2.all.sites.swe.mae.2018.png'),width=1000,height=700)
layout(mat = matrix(c(1,2), 
                    nrow = 2, 
                    ncol = 1),
             heights = c(1, 1.5),    # Heights of the two rows
             widths = 1)     # Widths of the two columns

##par(mfrow=c(2,1))

par(mar=c(0,5,3,3))
plot(0:alen,0:alen,xlab='',ylab='MAE (mm)',yaxs='i',
     col='white',main='',cex.axis=1.75,cex.lab=1.75,cex.main=2,
     xlim=c(1,alen),ylim=c(0,600),axes=FALSE)
##axis(1,at=1:alen,site.names[ranked.elevs],cex=1.75,cex.axis=1.75,las=2)
axis(2,at=c(0,200,400,600),c(0,200,400,600),cex=1.75,cex.axis=1.75)
abline(h=seq(0,1000,100),lty=2,col='gray',lwd=2)
abline(v=1:alen,col='gray')
for (j in 1:alen) {
    print(elevs[ranked.elevs[j]])
    boxplot(at=j-0.175,x=era.snodas.mae[ranked.elevs[j],],add=TRUE,axes=F,boxwex=0.7,col='blue',border='blue')
    boxplot(at=j+0.175,x=ncep2.snodas.mae[ranked.elevs[j],],add=TRUE,axes=F,boxwex=0.7,col='green',border='green')
}

text(x=16.6,y=550,'SNODAS',cex=1.25)
legend('topleft',leg=c('ERA-I','NCEP2'),col=c('blue','green'),pch=15,cex=1.5)

box(which='plot')

par(mar=c(10,5,0.1,3))
plot(0:alen,0:alen,xlab='',ylab='MAE (mm)',yaxs='i',
     col='white',main='',cex.axis=1.75,cex.lab=1.75,cex.main=2,
     xlim=c(1,alen),ylim=c(0,1200),axes=FALSE)
axis(1,at=1:alen,site.names[ranked.elevs],cex=1.75,cex.axis=1.75,las=2)
axis(2,at=seq(0,1200,200),seq(0,1200,200),cex=1.75,cex.axis=1.75)
abline(h=seq(0,2000,200),lty=2,col='gray',lwd=2)
abline(v=1:alen,col='gray')
for (j in 1:alen) {
    print(elevs[ranked.elevs[j]])
    boxplot(at=j-0.175,x=all.era.mae[ranked.elevs[j],],add=TRUE,axes=F,boxwex=0.7,col='blue',border='blue')
    boxplot(at=j+0.175,x=all.ncep2.mae[ranked.elevs[j],],add=TRUE,axes=F,boxwex=0.7,col='green',border='green')
    points(x=j,y=all.snodas.mae[ranked.elevs[j]],pch='-',cex=5,col='red')
}
text(x=16.0,y=1100,'Courses and Pillows',cex=1.25)
legend('topleft',leg=c('ERA-I','SNODAS','NCEP2'),col=c('blue','red','green'),pch=15,cex=1.5)

box(which='plot')

dev.off()
