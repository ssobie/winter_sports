##Script to compare snowdepth, SWE and snow cover between modelled values, 
##snow course sites and MODIS snow cover
library(ncdf4)
library(plotrix)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)
source('/storage/home/ssobie/code/repos/winter_sports/site.coordinates.r',chdir=T)

##---------------------------------------
##Observation data
read_course_obs <- function(site) {

   obs.file <- paste('/storage/data/projects/rci/data/winter_sports/obs/snow_courses/',site,'_snow_course.csv',sep='')
   obs.data <- read.csv(obs.file,header=T,as.is=T)
   obs.dates <- format(as.Date(obs.data[,1]),'%Y-%m-%d')
   obs.swe <- obs.data[,3] ##mm
   obs.na <- is.na(obs.swe)
   obs.swe <- obs.swe[!obs.na]
   obs.dates <- as.Date(obs.dates[!obs.na])

   rv <- list(dates=obs.dates,swe=obs.swe)
   return(rv)
}

read_snow_sim <- function(site,model,type) {

    model.dir <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow_series/'

    swe.file <- paste0(model.dir,site,'_',model,'_PNWNAmet_',type,'_snow_model_data.csv')
    swe.data <- read.csv(swe.file,header=T,as.is=T)
    swe.values <- swe.data$SWE*1000
    swe.dates <- as.Date(swe.data$Dates)
    rv <- list(dates=swe.dates,swe=swe.values)
    return(rv)
}


##-----------------------------------------

sites <- c('callaghan',
           'orchid_lake',
           'palisade_lake',
           'grouse_mountain',
           'dog_mountain',
           'stave_lake',
           'nahatlatch',
           'wahleach',
           'klesilkwa',
           'lightning_lake',
           'brookmere',
           'shovelnose_mountain',
           'hamilton_hill')


sites <- c('callaghan',
           'palisade_lake',
           'stave_lake',
           'nahatlatch',
           'klesilkwa',
           'brookmere')

site.names <- c('Callaghan','Palisade Lake',
                'Stave Lake','Nahatlatch',
                'Klesilkwa','Brookmere')


sites <- c('orchid_lake',           
           'grouse_mountain',
           'dog_mountain',                      
           'wahleach',           
           'lightning_lake',           
           'shovelnose_mountain',
           'hamilton_hill')

site.names <- c('Orchid Lake','Grouse Mountain',
                'Dog Mountain','Wahleach',
                'Lightning Lake','Shovelnose Mountain',
                'Hamilton Hill')



##Selected
sites <- c('grouse_mountain',
           'nahatlatch',
           'brookmere')

site.names <- c('Grouse Mountain',
           'Nahatlatch',
           'Brookmere')


##South Eastern Sites                
sites <- c('lightning_lake',
           'klesilkwa',
           'wahleach')

site.names <- c('Lightning Lake',
                'Klesilkwa',
                'Wahleach')


##Central Sites
sites <- c('dickson_lake',
           'stave_lake',
           'nahatlatch')

site.names <- c('Dickson Lake',
                'Stave Lake',
                'Nahatlatch')

##North Shore Sites
##           'grouse_mountain','Grouse Mountain',
sites <- c('orchid_lake',           
           'dog_mountain',
           'palisade_lake')                      
site.names <- c('Orchid Lake',
                'Dog Mountain','Palisade Lake')

##North Eastern Sites
sites <- c('brookmere',
           'shovelnose_mountain',
           'hamilton_hill')
           
site.names <- c('Brookmere',
                'Shovelnose Mountain',
                'Hamilton Hill')

##Northern Sites
sites <- c('duffey_lake',
           'mcgillivray_pass',
           'gnawed_mountain')

site.names <- c('Duffey Lake',
                'McGillivray Pass',
                'Gnawed Mountain')


type <- 'PRISM_TPS'

course.site.swe <- vector(mode='list',length=length(sites))
era5.site.swe <- vector(mode='list',length=length(sites))
pnw.site.swe <- vector(mode='list',length=length(sites))


##Loop over sites
model.dir <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow_series/'
plot.dir <- '/storage/data/projects/rci/data/winter_sports/plots/'

png(file=paste0(plot.dir,'pnw.era5.tps.northern.sites.swe.comparison.2020.png'),width=7,height=6,units='in',res=600,pointsize=6,bg='white')
par(mfrow=c(3,1),mar=c(0,0,0,0),oma=c(7,8,3,3))

for (i in seq_along(sites)) {
    site <- sites[i]
    print(site)

    coords <- get_coordinates(site)

    ##Snow Course Data
    course.file <- paste0('/storage/data/projects/rci/data/winter_sports/obs/snow_courses/',site,'_snow_course.csv',sep='')
    course.data <- read.csv(course.file,header=T,as.is=T)
    course.dates <- format(as.Date(course.data[,1]),'%Y-%m-%d')
    course.swe <- course.data[,3] ##mm
    course.pack <- course.data[,2] ##cm
    ##course.dense <-  course.data[,4]

    era5.swe.sim <- read_snow_sim(site,'ERA5',type)
    pnw.swe.sim <- read_snow_sim(site,'PNWNAmet',type)

    era5.course.subset <- format(as.Date(era5.swe.sim$dates),'%Y-%m-%d') %in% course.dates
    course.era5.subset <- course.dates %in% format(as.Date(era5.swe.sim$dates),'%Y-%m-%d')

    pnw.course.subset <- format(as.Date(pnw.swe.sim$dates),'%Y-%m-%d') %in% course.dates
    course.pnw.subset <- course.dates %in% format(as.Date(pnw.swe.sim$dates),'%Y-%m-%d')


    ##SNODAS Cell
    snodas.dir <- '/storage/data/projects/rci/data/winter_sports/obs/SNODAS/ncdf4_files/'
    snodas.file <- 'swe_snodas_modis_grid_van_whistler_20100101-20181231.nc'
    snc <- nc_open(paste0(snodas.dir,snodas.file))
    lon <- ncvar_get(snc,'lon')
    lat <- ncvar_get(snc,'lat')
    lon.ix <- which.min(abs(coords[1]-lon))
    lat.ix <- which.min(abs(coords[2]-lat))

    snodas.dates <- as.character(netcdf.calendar(snc))
    snodas.swe <- ncvar_get(snc,'swe',start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
    nc_close(snc)
 
##Snow Course Comparison
    yupp <- max(c(max(course.swe,na.rm=T),max(era5.swe.sim$swe),max(snodas.swe,na.rm=T),1000))
    ymax <- max(c(max(course.swe,na.rm=T)+100,max(pnw.swe.sim$swe),max(snodas.swe,na.rm=T)+100))


    plot(as.Date(course.dates),course.swe,cex=1.5,col='black',pch=16,yaxs='i',xaxs='i',
             xlim=c(as.Date('1979-08-01'),as.Date('2018-10-31')),ylim=c(0,ymax),
             main='',xlab='Date',ylab='SWE (mm)', cex.lab=1.95,axes=F)
    ##abline(h=c(500,1000,1500,2000,2500,3000),col='gray',lty=2,lwd=0.5)
    abline(h=seq(0,3600,300),col='gray',lty=2,lwd=0.5)

    if (i==3) {axis(1,at=as.Date(c('1980-01-01','1990-01-01','2000-01-01','2010-01-01','2018-01-01')),
           label=c('1980','1990','2000','2010','2018'),cex.axis=2.5,mgp=c(3,2,0))}
    axis(2,at=seq(0,yupp,round((yupp-100)/3,-2)),label=seq(0,yupp,round((yupp-100)/3,-2)),cex.axis=2.5,mgp=c(3,2,0))

    lines(as.Date(era5.swe.sim$dates),era5.swe.sim$swe,lwd=3,col='blue')
    lines(as.Date(pnw.swe.sim$dates),pnw.swe.sim$swe,lwd=2,col='green')
    lines(as.Date(snodas.dates),snodas.swe,lwd=1.5,col='red')
    points(as.Date(course.dates),course.swe,cex=1.75,col='black',pch=16)
    text(as.Date('2010-01-01'),0.95*ymax,site.names[i],cex=2.5)
    if (i==3) {
       legend('topleft',legend=c('Course Obs.','SNODAS','ERA5','PNWNAmet'),col=c('black','red','blue','green'),pch=15,cex=1.75,pt.cex=3)
    }
    box(which='plot')
    abline(h=0)

}        

mtext("Date",side=1,outer=TRUE,cex=2.0,line=4.6)
mtext("Snow Water Equivalent (mm)",side=2,outer=TRUE,cex=2.0,line=4.6)


dev.off()

