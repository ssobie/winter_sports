##Script to compare snowdepth, SWE and snow cover between modelled values, 
##snow course sites and MODIS snow cover
library(ncdf4)
library(plotrix)

##source('/storage/home/ssobie/code/repos/winter_sports/test.snow.model.r',chdir=T)
source('/storage/home/ssobie/code/repos/winter_sports/site.coordinates.r',chdir=T)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)


get_mae_values <- function(sites,asps) {

course.site.swe <- vector(mode='list',length=length(sites))
pillow.site.swe <- vector(mode='list',length=length(asps))
model.site.swe <- vector(mode='list',length=length(sites))

era5.file <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/grouse_mountain_ERA5_800m_data.csv'
era5.data <- read.csv(era5.file,header=T,as.is=T)
era5.dates <- era5.data$Dates

pnwnamet.file <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/grouse_mountain_PNWNAmet_800m_data.csv'
pnwnamet.data <- read.csv(pnwnamet.file,header=T,as.is=T)
pnwnamet.dates <- pnwnamet.data$Dates

pnwnamet.course.bias <- era5.course.bias <- rep(0,length(sites))
pnwnamet.pillow.bias <- era5.pillow.bias <- rep(0,length(asps))
pnwnamet.snodas.bias <- era5.snodas.bias <- rep(0,length(sites)+length(asps))
snodas.site.bias <- rep(0,length(sites))
snodas.asp.bias <- rep(0,length(asps))

model.dir <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow_series/'

snodas.file <- "/storage/data/projects/rci/data/winter_sports/obs/SNODAS/ncdf4_files/swe_snodas_prism_grid_van_whistler_20100101-20181231.nc"
snc <- nc_open(snodas.file)
lon <- ncvar_get(snc,'lon')
lat <- ncvar_get(snc,'lat')
snodas.dates <- as.character(netcdf.calendar(snc))

pnwnamet.snodas.match <-  format(as.Date(pnwnamet.dates),'%Y-%m-%d') %in% format(as.Date(snodas.dates),'%Y-%m-%d') 
snodas.pnwnamet.match <-  format(as.Date(snodas.dates),'%Y-%m-%d') %in% format(as.Date(pnwnamet.dates),'%Y-%m-%d') 

era5.snodas.match <-  format(as.Date(era5.dates),'%Y-%m-%d') %in% format(as.Date(snodas.dates),'%Y-%m-%d') 
snodas.era5.match <-  format(as.Date(snodas.dates),'%Y-%m-%d') %in% format(as.Date(era5.dates),'%Y-%m-%d') 

for (i in seq_along(sites)) {
    site <- sites[i]
    print(site)
    coords <- get_coordinates(site)

    ##Snow Course Data
    course.file <- paste('/storage/data/projects/rci/data/winter_sports/obs/snow_courses/',site,'_snow_course.csv',sep='')
    course.data <- read.csv(course.file,header=T,as.is=T)
    course.dates <- format(as.Date(course.data[,1]),'%Y-%m-%d')
    course.swe <- course.data[,3] ##mm
    era5.course.subset <- format(as.Date(era5.dates),'%Y-%m-%d') %in% course.dates
    course.era5.subset <- course.dates %in% format(as.Date(era5.dates),'%Y-%m-%d')
    pnwnamet.course.subset <- format(as.Date(pnwnamet.dates),'%Y-%m-%d') %in% course.dates
    course.pnwnamet.subset <- course.dates %in% format(as.Date(pnwnamet.dates),'%Y-%m-%d')

    ##SNODAS Data at Courses
    lon.ix <- which.min(abs(coords[1]-lon))
    lat.ix <- which.min(abs(coords[2]-lat))
    snodas.swe <- ncvar_get(snc,'swe',start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
    snodas.date.subset <- format(as.Date(snodas.dates),'%Y-%m-%d') %in% course.dates
    snodas.course.subset <- course.dates %in% format(as.Date(snodas.dates),'%Y-%m-%d')
    snodas.site.bias[i] <- mean((snodas.swe[snodas.date.subset]-course.swe[snodas.course.subset]),na.rm=T)
      
    swe.era5 <- read.csv(paste0(model.dir,site,'_ERA5_PNWNAmet_PRISM_TPS_snow_model_data.csv'),header=T,as.is=T)
    era5.swe.mean <- swe.era5$SWE*1000 ###apply(swe.era5,1,mean,na.rm=T)
    swe.pnwnamet <- read.csv(paste0(model.dir,site,'_PNWNAmet_PNWNAmet_PRISM_TPS_snow_model_data.csv'),header=T,as.is=T)
    pnwnamet.swe.mean <- swe.pnwnamet$SWE*1000 ##apply(swe.pnwnamet,1,mean,na.rm=T)

    era5.course.bias[i] <- mean((era5.swe.mean[era5.course.subset]-course.swe[course.era5.subset]),na.rm=T)
    pnwnamet.course.bias[i] <- mean((pnwnamet.swe.mean[pnwnamet.course.subset]-course.swe[course.pnwnamet.subset]),na.rm=T)
    era5.snodas.bias[i] <- mean((era5.swe.mean[era5.snodas.match]-snodas.swe[snodas.era5.match]),na.rm=T)      
    pnwnamet.snodas.bias[i] <- mean((pnwnamet.swe.mean[pnwnamet.snodas.match]-snodas.swe[snodas.pnwnamet.match]),na.rm=T)      
    ##if(site=='palisade_lake'){browser()}
}        

for (i in seq_along(asps)) {
    asp <- asps[i]
    coords <- get_coordinates(asp)

    ##Snow Pillow Data
    pillow.file <- paste('/storage/data/projects/rci/data/winter_sports/obs/snow_pillow/',asp,'.csv',sep='')
    pillow.data <- read.csv(pillow.file,header=T,as.is=T)
    pillow.dates <- format(as.Date(pillow.data[,2]),'%Y-%m-%d')
    pillow.swe <- pillow.data[,11] ##mm
    era5.pillow.subset <- format(as.Date(era5.dates),'%Y-%m-%d') %in% pillow.dates
    pillow.era5.subset <- pillow.dates %in% format(as.Date(era5.dates),'%Y-%m-%d')
    pnwnamet.pillow.subset <- format(as.Date(pnwnamet.dates),'%Y-%m-%d') %in% pillow.dates
    pillow.pnwnamet.subset <- pillow.dates %in% format(as.Date(pnwnamet.dates),'%Y-%m-%d')

    ##SNODAS Data at Courses
    lon.ix <- which.min(abs(coords[1]-lon))
    lat.ix <- which.min(abs(coords[2]-lat))
    snodas.swe <- ncvar_get(snc,'swe',start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
    snodas.date.subset <- format(as.Date(snodas.dates),'%Y-%m-%d') %in% pillow.dates
    snodas.pillow.subset <- pillow.dates %in% format(as.Date(snodas.dates),'%Y-%m-%d')

    swe.era5 <- read.csv(paste0(model.dir,asp,'_ERA5_PRISM_TPS_snow_model_data.csv'),header=T,as.is=T)
    era5.swe.mean <- swe.era5$SWE*1000 ##apply(swe.era5,1,mean,na.rm=T)
    swe.pnwnamet <- read.csv(paste0(model.dir,asp,'_PNWNAmet_PRISM_TPS_snow_model_data.csv'),header=T,as.is=T)
    pnwnamet.swe.mean <- swe.pnwnamet$SWE*1000 ##apply(swe.pnwnamet,1,mean,na.rm=T)

    snodas.asp.bias[i] <- mean((snodas.swe[snodas.date.subset]-pillow.swe[snodas.pillow.subset]),na.rm=T)      
    era5.pillow.bias[i] <- mean((era5.swe.mean[era5.pillow.subset]-pillow.swe[pillow.era5.subset]),na.rm=T)
    pnwnamet.pillow.bias[i] <- mean((pnwnamet.swe.mean[pnwnamet.pillow.subset]-pillow.swe[pillow.pnwnamet.subset]),na.rm=T)
    era5.snodas.bias[i+ilen] <- mean((era5.swe.mean[era5.snodas.match]-snodas.swe[snodas.era5.match]),na.rm=T)      
    pnwnamet.snodas.bias[i+ilen] <- mean((pnwnamet.swe.mean[pnwnamet.snodas.match]-snodas.swe[snodas.pnwnamet.match]),na.rm=T)      


}        

nc_close(snc)

all.era5.bias <- c(era5.course.bias,era5.pillow.bias)
all.pnwnamet.bias <- c(pnwnamet.course.bias,pnwnamet.pillow.bias)
all.snodas.bias <- c(snodas.site.bias,snodas.asp.bias)

rv <- list(era=all.era5.bias,pnw=all.pnwnamet.bias,sno=all.snodas.bias)
           
return(rv)

}

##*******************************************************************************


sites <- c('shovelnose_mountain','brookmere','lightning_lake','callaghan','orchid_lake',
           'palisade_lake','grouse_mountain','dog_mountain','stave_lake','nahatlatch',
           'wahleach','klesilkwa','hamilton_hill','dickson_lake','disappointment_lake',
           'duffey_lake','gnawed_mountain','highland_valley','mcgillivray_pass',
           'sumallo_river_west','great_bear')
ilen <- length(sites)
asps <- c('upper_squamish','spuzzum_creek','chilliwack_river','tenquille_lake',
           'wahleach_lake','blackwall_peak_pillow')

site.names <- c('Shovelnose', ##\nMountain',
                'Brookmere',
                'Lightning', ##\nLake',
                'Callaghan',
                'Orchid L.', ##\nLake',
                'Palisade', ##\nLake',
                'Grouse', ##\nMountain',
                'Dog Mtn.', ##\nMountain',
                'Stave L.', ##\nLake',
                'Nahatlatch',
                'Wahleach',
                'Klesilkwa',
                'Hamilton', ##\nHill',
                'Dickson', ##\nLake',
                'Disappoint.', ##ment', ##\nLake',
                'Duffey L.', ##\nLake',
                'Gnawed', ##\nMountain',
                'Highland', ##\nValley',
                'Mcgillivray', ##\nPass',
                'Sumallo W.', ##\nWest',          
                'Bear', ##(Great)
                'Squamish', ##Upper\n
                'Spuzzum', ##\nCreek',
                'Chilliwack', ##\nRiver',
                'Tenquille L.', ##\nLake',
                'Wahleach L.', ##\nLake',
                'Blackwall P.') ##\nPeak')

course.lons <- rep(0,length(sites))
pillow.lons <- rep(0,length(asps))
course.lats <- rep(0,length(sites))
pillow.lats <- rep(0,length(asps))

course.elevs <- rep(0,length(sites))
pillow.elevs <- rep(0,length(asps))

for (i in seq_along(sites)) {
    site <- sites[i]
    coords <- get_coordinates(site)
    course.lons[i] <- coords[1]
    course.lats[i] <- coords[2]
    course.elevs[i] <- coords[3]
}

for (i in seq_along(asps)) {
    asp <- asps[i]
    coords <- get_coordinates(asp)
    pillow.lons[i] <- coords[1]
    pillow.lats[i] <- coords[2]
    pillow.elevs[i] <- coords[3]
}

lons <- c(course.lons,pillow.lons)
lats <- c(course.lats,pillow.lats)
elevs <- c(course.elevs,pillow.elevs)
alen <- length(sites)+length(asps)

cal.mae <- get_mae_values(sites,asps) 


##-------------------------------------------------------------------
if (1==0) {
snodas.bias <- list(lons=lons,lats=lats,bias=cal.mae$sno)

##save(snodas.bias,file=paste0('/storage/data/projects/rci/data/winter_sports/snodas.bias.RData'))

##browser()


##Loop over sites
plot.dir <- '/storage/data/projects/rci/data/winter_sports/plots/'

leg.title <- 'mm'
ranked.lons <- order(lons)
ranked.elevs <- order(elevs)

##png(filename=paste0(plot.dir,'pnwnamet.era5.tps.cal.sites.swe.bias.2020.png'),width=1200,height=700)
png(file=paste0(plot.dir,'pnwnamet.era5.tps.cal.sites.swe.bias.45.png'),width=10,height=5,units='in',res=600,pointsize=6,bg='white')

layout(mat = matrix(c(1,2),
                    nrow = 2,
                    ncol = 1),
             heights = c(1.6, 1.4),    # Heights of the two rows
             widths = 1)     # Widths of the two columns

par(mar=c(7,5,2,3))
plot(0:alen,0:alen,xlab='',ylab='SWE Bias (mm)',yaxs='i',
     col='white',main='',cex.axis=1.75,cex.lab=1.75,cex.main=2,
     xlim=c(1,alen),ylim=c(-1200,800),axes=FALSE)
###axis(1,at=1:alen,site.names[ranked.lons],cex=1.75,cex.axis=1.75,las=2)
axis(1,at=1:alen,labels=FALSE,cex=1.75,cex.axis=1.75,las=2)
par(xpd=NA)
text((1:alen)-0.8,par("usr")[3] - 260, labels=site.names[ranked.lons],srt=45,pos=1,cex=1.75)
par(xpd=FALSE)
print(par("usr"))

axis(2,at=c(-1200,-400,400,800),label=c(-1200,-400,400,800),cex=1.75,cex.axis=1.75)
abline(h=seq(-800,400,400),lty=2,col='gray',lwd=2)
abline(v=1:alen,col='gray')
abline(h=0,lwd=1.5)
for (j in 1:alen) {
    print(lons[ranked.lons[j]])
    ##boxplot(at=j-0.175,x=cal.mae$era[ranked.lons[j]],add=TRUE,axes=F,boxwex=0.7,col='blue',border='blue')
    points(x=j-0.175,y=cal.mae$era[ranked.lons[j]],col='blue',pch=18,cex=2.25)
    points(x=j+0.175,y=cal.mae$pnw[ranked.lons[j]],col='green',pch=18,cex=2.25)
    points(x=j,y=cal.mae$sno[ranked.lons[j]],pch='-',cex=5,col='red')
}
##text(x=16.6,y=550,'SNODAS',cex=1.25)

legend('bottomright',leg=c('ERA5','PNWNAmet','SNODAS'),col=c('blue','green','red'),pch=c(18,18,15),
       pt.cex=c(2,2,1.5),box.lwd=1.5,cex=1.5)

box(which='plot',lwd=1.5)

}

sites <- c('blackwall_peak_course','boston_bar_lower','boston_bar_upper','burwell_lake',
           'chapman_creek','cornwall_hills','diamond_head','edwards_lake',
           'hollyburn','hope',
           'loch_lomond','lytton','mount_seymour','new_tashme',
           'ottomite','pavilion_mountain',
           'sumallo_river','tenquille_course','whistler_mountain','wolverine_creek' )

elen <- length(sites)

site.names <- c('Blackwall',
                'Boston U.',
                'Boston L.',
                'Burwell',
                'Chapman', 
                'Cornwall',
                'Diamond', 
                'Edwards', 
                'Hollyburn',
                'Hope',
                'Lomond',
                'Lytton', 
                'Seymour',
                'Tashme',
                'Ottomite',
                'Pavilion',
                'Sumallo',
                'Tenquille',
                'Whistler', 
                'Wolverine')

course.lons <- rep(0,length(sites))
course.lats <- rep(0,length(sites))

course.elevs <- rep(0,length(sites))

for (i in seq_along(sites)) {
    site <- sites[i]
    coords <- get_coordinates(site)
    course.lons[i] <- coords[1]
    course.lats[i] <- coords[2]
    course.elevs[i] <- coords[3]
}

lons <- course.lons
lats <- course.lats
elevs <- course.elevs
elen <- length(sites)
ranked.lons <- order(lons)
ranked.elevs <- order(elevs)

val.mae <- get_mae_values(sites,asps) 
slen <- length(sites)
val.mae <- lapply(val.mae,function(x,slen){x[1:slen]},slen)

browser()

par(mar=c(7,5,2,9))
plot(0:elen,0:elen,xlab='',ylab='SWE Bias (mm)',yaxs='i',
     col='white',main='',cex.axis=1.75,cex.lab=1.75,cex.main=2,
     xlim=c(1,elen),ylim=c(-800,800),axes=FALSE)
###axis(1,at=1:elen,site.names[ranked.lons],cex=1.75,cex.axis=1.75,las=2)
###axis(1,at=1:elen,site.names[ranked.lons],cex=1.75,cex.axis=1.75,las=2)
axis(1,at=1:elen,labels=FALSE,cex=1.75,cex.axis=1.75,las=2)
par(xpd=NA)
text((1:elen)-0.7,par("usr")[3] - 260, labels=site.names[ranked.lons],srt=45,pos=1,cex=1.75)
par(xpd=FALSE)

axis(2,at=c(-800,-400,400,800),label=c(-800,-400,400,800),cex=1.75,cex.axis=1.75)
abline(h=seq(-400,400,400),lty=2,col='gray',lwd=2)
abline(v=1:alen,col='gray')
abline(h=0,lwd=1.5)
for (j in 1:elen) {
    print(lons[ranked.lons[j]])
    ##boxplot(at=j-0.175,x=val.mae$era[ranked.lons[j]],add=TRUE,axes=F,boxwex=0.7,col='blue',border='blue')
    points(x=j-0.175,y=val.mae$era[ranked.lons[j]],col='blue',pch=18,cex=2.25)
    points(x=j+0.175,y=val.mae$pnw[ranked.lons[j]],col='green',pch=18,cex=2.25)
    ##points(x=j,y=val.mae$sno[ranked.lons[j]],pch='-',cex=5,col='red')
}
legend('bottomright',leg=c('ERA5','PNWNAmet'),col=c('blue','green'),pch=c(18,18),
       pt.cex=c(2,2,1.5),box.lwd=1.5,cex=1.5)

##par(xpd=NA)
##legend('topright',inset=c(-0.08,0),leg=c('ERA5','PNW','SNO'),col=c('blue','green','red'),pch=c(18,18,15),pt.cex=c(2,2,1.5),cex=1.5,bo##x.lwd=1.5)
##par(xpd=FALSE)

box(which='plot',lwd=1.5)


dev.off()

browser()

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
##axis(1,at=1:alen,site.names[ranked.lons],cex=1.75,cex.axis=1.75,las=2)
axis(2,at=c(0,200,400,600),c(0,200,400,600),cex=1.75,cex.axis=1.75)
abline(h=seq(0,1000,100),lty=2,col='gray',lwd=2)
abline(v=1:alen,col='gray')
for (j in 1:alen) {
    print(lons[ranked.lons[j]])
    boxplot(at=j-0.175,x=era.snodas.mae[ranked.lons[j],],add=TRUE,axes=F,boxwex=0.7,col='blue',border='blue')
    boxplot(at=j+0.175,x=ncep2.snodas.mae[ranked.lons[j],],add=TRUE,axes=F,boxwex=0.7,col='green',border='green')
}

text(x=16.6,y=550,'SNODAS',cex=1.25)
legend('topleft',leg=c('ERA-I','NCEP2'),col=c('blue','green'),pch=15,cex=1.5)

box(which='plot')

par(mar=c(10,5,0.1,3))
plot(0:alen,0:alen,xlab='',ylab='MAE (mm)',yaxs='i',
     col='white',main='',cex.axis=1.75,cex.lab=1.75,cex.main=2,
     xlim=c(1,alen),ylim=c(0,1200),axes=FALSE)
axis(1,at=1:alen,site.names[ranked.lons],cex=1.75,cex.axis=1.75,las=2)
axis(2,at=seq(0,1200,200),seq(0,1200,200),cex=1.75,cex.axis=1.75)
abline(h=seq(0,2000,200),lty=2,col='gray',lwd=2)
abline(v=1:alen,col='gray')
for (j in 1:alen) {
    print(lons[ranked.lons[j]])
    boxplot(at=j-0.175,x=all.era.mae[ranked.lons[j],],add=TRUE,axes=F,boxwex=0.7,col='blue',border='blue')
    boxplot(at=j+0.175,x=all.ncep2.mae[ranked.lons[j],],add=TRUE,axes=F,boxwex=0.7,col='green',border='green')
    points(x=j,y=all.snodas.mae[ranked.lons[j]],pch='-',cex=5,col='red')
}
text(x=16.0,y=1100,'Courses and Pillows',cex=1.25)
legend('topleft',leg=c('ERA-I','SNODAS','NCEP2'),col=c('blue','red','green'),pch=15,cex=1.5)

box(which='plot')

dev.off()
