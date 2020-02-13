##Script to compare snowdepth, SWE and snow cover between modelled values, 
##snow course sites and MODIS snow cover
library(ncdf4)
library(plotrix)
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

##North Shore Sites
sites <- c('orchid_lake',           
           'grouse_mountain',
           'dog_mountain',
           'palisade_lake')                      
site.names <- c('Orchid Lake','Grouse Mountain',
                'Dog Mountain','Palisade Lake')


course.site.swe <- vector(mode='list',length=length(sites))
model.site.swe <- vector(mode='list',length=length(sites))


##Loop over sites

plot.dir <- '/storage/data/projects/rci/data/winter_sports/ncc_2019/plots/'
png(file=paste0(plot.dir,'north.shore.sites.swe.courses.only.2019.png'),width=5,height=4,units='in',res=600,pointsize=6,bg='white')
par(mfrow=c(4,1))    
par(mar=c(0,6.1,0,0),oma=c(6,0,4,4))
par(mgp=c(4,1.5,0))

x.dates <- as.Date(c('1980-01-01','1985-01-01','1990-01-01','1995-01-01',
                     '2000-01-01','2005-01-01','2010-01-01','2015-01-01'))
x.names <- as.character(format(x.dates,'%Y'))

for (i in seq_along(sites)) {
    site <- sites[i]
    print(site)

    ##Snow Course Data
    course.file <- paste0('/storage/data/projects/rci/data/winter_sports/obs/snow_courses/',site,'_snow_course.csv',sep='')
    ##course.file <- paste('/storage/data/projects/rci/data/assessments/snow_model)
    course.data <- read.csv(course.file,header=T,as.is=T)
    course.dates <- format(as.Date(course.data[,1]),'%Y-%m-%d')
    course.swe <- course.data[,3] ##mm
    course.pack <- course.data[,2] ##cm
    course.dense <-  course.data[,4]

##Snow Course Comparison
    yupp <- max(c(max(course.swe,na.rm=T),1000))
    ymax <- max(course.swe*1.2,na.rm=T)
    ##par(mar=c(5.1,5,2.1,2.1))
    plot(as.Date(course.dates),course.swe,cex=1.1,col='blue',pch='*',
             xlim=c(as.Date('1981-01-01'),as.Date('2017-04-30')),ylim=c(0,ymax),yaxs='i',
             main='',xlab='Date',ylab='SWE (mm)', cex.lab=1.95,axes=F)
    axis(2,at=seq(0,yupp,round((yupp-100)/3,-2)),label=seq(0,yupp,round((yupp-100)/3,-2)),cex.axis=1.45)
    abline(v=x.dates,col='gray',lwd=0.5,lty=2)
    points(as.Date(course.dates),course.swe,cex=1.7,col='blue',pch=16)

    text(as.Date('1986-01-01'),0.9*ymax,site.names[i],cex=2)
    
    box(which='plot')
    abline(h=0)
}        
  axis(1,at=x.dates,label=x.names,cex.axis=1.95)
  mtext("Date",side=1,outer=TRUE,cex=1.5,line=3.6)

dev.off()

