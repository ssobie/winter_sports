##Script to compare snowdepth, SWE and snow cover between modelled values, 
##snow course sites and MODIS snow cover
library(ncdf4)
library(plotrix)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##source('/storage/home/ssobie/code/repos/winter_sports/test.snow.model.r',chdir=T)

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

hyper.snow <- function(pr,tasmax,tasmin,coeffs) {
        
        tas <- (tasmax+tasmin)/2                # degrees C
        frac <- coeffs$a*(tanh(coeffs$b*(tas-coeffs$c))-coeffs$d)
        sample <- runif(length(tas),min=0,max=100)
        test <- sample > frac        
        snow.type <- rep(TRUE,length(tas))
        snow.type[test] <- FALSE

        NewSnowWatEq <- pr*0.001
        NewSnowWatEq[!snow.type] <- 0
        R_m <- pr*0.001
        R_m[snow.type] <- 0
        rv <- list(snow=NewSnowWatEq,
                   rain=R_m)
        return(rv)
}


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



##North Eastern Sites
sites <- c('brookmere',
           'shovelnose_mountain',
           'hamilton_hill')
           
site.names <- c('Brookmere',
                'Shovelnose Mountain',
                'Hamilton Hill')

##South Eastern Sites                
sites <- c('lightning_lake',
           'klesilkwa',
           'wahleach')

site.names <- c('Lightning Lake',
                'Klesilkwa',
                'Wahleach')

##Central Sites
sites <- c('callaghan',
           'stave_lake',
           'nahatlatch')

site.names <- c('Callaghan',
                'Stave Lake',
                'Nahatlatch')

##North Shore Sites
sites <- c('orchid_lake',           
           'grouse_mountain',
           'dog_mountain',
           'palisade_lake')                      
site.names <- c('Orchid Lake','Grouse Mountain',
                'Dog Mountain','Palisade Lake')

##Selected
sites <- c('grouse_mountain',
           'nahatlatch',
           'brookmere')

site.names <- c('Grouse Mountain',
           'Nahatlatch',
           'Brookmere')


model <- 'ERA'

course.site.swe <- vector(mode='list',length=length(sites))
model.site.swe <- vector(mode='list',length=length(sites))

slen <- 1001

model.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/grouse_mountain_',model,'_800m_data.csv')
model.data <- read.csv(model.file,header=T,as.is=T)

swe.sims <- matrix(0,nrow=slen,ncol=dim(model.data)[1])
snow.sims <- matrix(0,nrow=slen,ncol=dim(model.data)[1])


##Loop over sites
model.dir <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sims/'
plot.dir <- '/storage/data/projects/rci/data/winter_sports/plots/course_comparison/'
##png(file=paste0(plot.dir,model,'.subset.sites.swe.comparison.2018.png'),width=1400,height=1000)
png(file=paste0(plot.dir,model,'.selected.sites.swe.comparison.2018.png'),width=1000,height=900)
par(mfrow=c(3,1))    

for (i in seq_along(sites)) {
    site <- sites[i]
    print(site)
    ##Reanalysis 800m data
    clim.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/',site,'_',model,'_800m_data.csv')
    clim.data <- read.csv(clim.file,header=T,as.is=T)
    pr.data <- clim.data$Pr
    tasmax.data <- clim.data$Tasmax
    tasmin.data <- clim.data$Tasmin
    tas.data <- clim.data$Tas
    dates <- clim.data$Dates

    coords <- get.coordinates(site)
    lat.bnds <- coords[2]
    elev <- coords[3]
    site.slope <- bc.slopes[which.min(abs(coords[1]-bc.lon)),which.min(abs(coords[2]-bc.lat))]    
    site.aspect <- bc.aspects[which.min(abs(coords[1]-bc.lon)),which.min(abs(coords[2]-bc.lat))]    

    ##Snow Course Data
    course.file <- paste0('/storage/data/projects/rci/data/winter_sports/obs/snow_courses/',site,'_snow_course.csv',sep='')
    ##course.file <- paste('/storage/data/projects/rci/data/assessments/snow_model)
    course.data <- read.csv(course.file,header=T,as.is=T)
    course.dates <- format(as.Date(course.data[,1]),'%Y-%m-%d')
    course.swe <- course.data[,3] ##mm
    course.pack <- course.data[,2] ##cm
    course.dense <-  course.data[,4]
    sb <- 1:length(course.dates) ##190:210
    date.subset <- format(as.Date(dates),'%Y-%m-%d') %in% course.dates
    course.subset <- course.dates %in% format(as.Date(dates),'%Y-%m-%d')

    swe.sims <- read.csv(paste0(model.dir,site,'.',tolower(model),'.swe.1001.csv'),header=T,as.is=T)
    swe.mean <- apply(swe.sims,1,mean,na.rm=T)

    ##SNODAS Cell
    snodas.dir <- '/storage/data/projects/rci/data/winter_sports/obs/SNODAS/ncdf4_files/'
    snodas.file <- 'swe_snodas_modis_grid_van_whistler_20100101-20181231.nc'
    snc <- nc_open(paste0(snodas.dir,snodas.file))
    lon <- ncvar_get(snc,'lon')
    lat <- ncvar_get(snc,'lat')
    lon.ix <- which.min(abs(coords[1]-lon))
    lat.ix <- which.min(abs(coords[2]-lat))
    if (grepl('(palisade_lake|wahleach)',site))
       lat.ix <- lat.ix+1

    snodas.dates <- as.character(netcdf.calendar(snc))
    snodas.swe <- ncvar_get(snc,'swe',start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
    nc_close(snc)

 
    ###all.ae <- rep(0,slen)
    ###for (j in 1:slen) {
    ###  all.ae[j] <- mean(abs(course.swe[course.subset]-swe.sims[date.subset,j]),na.rm=T)
    ###}  
    ##print(range(all.ae))
    ##print(mean(all.ae))

##Snow Course Comparison
    yupp <- max(c(max(course.swe,na.rm=T),max(swe.sims),max(snodas.swe,na.rm=T),1000))
    ymax <- max(c(max(course.swe,na.rm=T),max(swe.sims),max(snodas.swe,na.rm=T)))

    par(mar=c(5.1,5,2.1,2.1))
    plot(as.Date(course.dates)[sb],course.swe[sb],cex=1.5,col='black',pch=16,
             xlim=c(as.Date('1981-01-01'),as.Date('2017-04-30')),ylim=c(0,ymax),
             main='',xlab='Date',ylab='SWE (mm)', cex.lab=1.95,axes=F)
    axis(1,at=as.Date(c('1980-01-01','1990-01-01','2000-01-01','2010-01-01')),label=c('1980','1990','2000','2010'),cex.axis=1.95)
    axis(2,at=seq(0,yupp,round((yupp-100)/3,-2)),label=seq(0,yupp,round((yupp-100)/3,-2)),cex.axis=1.95)
    apply(swe.sims,2,function(x,y){lines(y,x,col='lightblue',lwd=2)},as.Date(dates))
    lines(as.Date(dates),swe.mean,lwd=3,col='blue')
    lines(as.Date(snodas.dates),snodas.swe,lwd=2,col='red')
    points(as.Date(course.dates),course.swe,cex=1.75,col='black',pch=16)
    text(as.Date('2015-01-01'),0.95*ymax,site.names[i],cex=2)
    if (i==3) {
       legend('topleft',legend=c('Course Obs.','SNODAS','Model','Model Mean'),col=c('black','red','lightblue','blue'),pch=16,cex=1.75)
    }
    box(which='plot')
    abline(h=0)


}        

dev.off()

