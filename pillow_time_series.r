##Script to compare snowdepth, SWE and snow cover between modelled values, 
##snow course sites and MODIS snow cover
library(ncdf4)
library(plotrix)
library(scales)

source('/storage/home/ssobie/code/repos/winter_sports/site.coordinates.r',chdir=T)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##Slope and Aspect Values
model <- 'PNWNAmet'
reanalysis <- 'PNWNAmet'
var.name <- 'SWE'
type <- 'PRISM_TPS'

sites <- c('blackwall_peak_pillow','spuzzum_creek','upper_squamish','chilliwack_river','tenquille_lake','wahleach_lake')
site.names <- c('Blackwall Peak','Spuzzum Creek','Upper Squamish','Chilliwack River','Tenquille Lake','Wahleach Lake')

    plot.dir <- '/storage/data/projects/rci/data/winter_sports/plots/'
    ##png(file=paste0(plot.dir,model,'.SWE.normalized.pillow.series.1001.2018.png'),width=1000,height=900)
    png(file=paste0(plot.dir,'PNWNAmet.ERA5.TPS.Elevation.SWE.normalized.pillow.series.reordered.png'),width=6,height=6,units='in',res=600,pointsize=6,bg='white')
    par(mfrow=c(6,1),mar=c(0,0,0,0),oma=c(7,8,3,3))    

##Loop over sites

for (i in seq_along(sites)) {
    site <- sites[i]
    print(site)
    clim.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/',site,'_ERA5_800m_data.csv')
    clim.data <- read.csv(clim.file,header=T,as.is=T)
    era5.dates <- clim.data$Dates

    clim.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/',site,'_PNWNAmet_800m_data.csv')
    clim.data <- read.csv(clim.file,header=T,as.is=T)
    pnw.dates <- clim.data$Dates

    snow.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow_series/',site,'_ERA5_PNWNAmet_',type,'_snow_model_data.csv')
    snow.sims <- read.csv(snow.file,header=T,as.is=T)
    era5.swe.sims <- snow.sims$SWE*1000

    snow.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow_series/',site,'_PNWNAmet_PNWNAmet_',type,'_snow_model_data.csv')
    snow.sims <- read.csv(snow.file,header=T,as.is=T)
    pnw.swe.sims <- snow.sims$SWE*1000

    coords <- get_coordinates(site)
    lat.bnds <- coords[2]
    elev <- coords[3]

    ##Snow Pillow Data
    pillow.file <- paste('/storage/data/projects/rci/data/winter_sports/obs/snow_pillow/',site,'.csv',sep='')
    pillow.data <- read.csv(pillow.file,header=T,as.is=T)
    pillow.dates <- format(as.Date(pillow.data[,2]),'%Y-%m-%d')
    pillow.tasmax <- pillow.data[,3]
    pillow.tasmin <- pillow.data[,5]
    pillow.tas <- (pillow.tasmax + pillow.tasmin)/2
    pillow.precip <- pillow.data[,7]##mm
    pillow.swe <- pillow.data[,11] ##mm
    pillow.pack <- pillow.data[,13] ##cm


    ##SNODAS Cell
    snodas.dir <- '/storage/data/projects/rci/data/winter_sports/obs/SNODAS/ncdf4_files/'    
    snodas.file <- 'swe_snodas_prism_grid_van_whistler_20100101-20181231.nc'
    snc <- nc_open(paste0(snodas.dir,snodas.file))
    lon <- ncvar_get(snc,'lon')
    lat <- ncvar_get(snc,'lat')    
    lon.ix <- which.min(abs(coords[1]-lon))
    lat.ix <- which.min(abs(coords[2]-lat))

    snodas.dates <- as.character(netcdf.calendar(snc))

    snodas.swe <- ncvar_get(snc,'swe',start=c(lon.ix,lat.ix,1),count=c(1,1,-1))
    nc_close(snc)    

    ##Normalized Series

      
    if (1==1) {
    ##SWE
    if (var.name=='SWE') {
      plot(as.Date(era5.dates),era5.swe.sims,xlab='Date',ylab='SWE (mm)',yaxs='i',
           type='l',lwd=3,col='blue',main='',cex.axis=1.75,cex.lab=1.75,cex.main=2,xaxs='i',axes=F,
           xlim=c(as.Date('1989-08-01'),as.Date('2018-10-31')),ylim=c(0,3100))
      abline(h=c(500,1000,1500,2000,2500,3000),col='gray',lty=2,lwd=0.5)

      axis(2,at=c(0,500,1000,1500,2000,2500,3000),labels=c('','','','','','',''),cex=2,mgp=c(3,2,0))
      axis(2,at=c(0,1500),label=c('0','1500'),cex.axis=2.5,mgp=c(3,2,0))
      if (i==1) {axis(2,at=c(0,1500,3000),label=c('0','1500','3000'),cex.axis=2.5,mgp=c(3,2,0))}
      points(as.Date(snodas.dates),snodas.swe,col=alpha('red',0.5),cex=1.1)
      points(as.Date(pillow.dates),pillow.swe,pch=16,col='black',cex=1)
      lines(as.Date(era5.dates),era5.swe.sims,col='blue',lwd=2.35)
      lines(as.Date(pnw.dates),pnw.swe.sims,col='green',lwd=1.85)
      text(as.Date('1995-01-01'),2600,site.names[i],cex=2.5)
      if (i==3) {legend('topright',legend=c('ASP Obs.'),col=c('black'),pch=15,cex=1.75,pt.cex=3)}
      if (i==4) {legend('topright',legend=c('SNODAS'),col=c('red'),pch=15,cex=1.75,pt.cex=3)}
      if (i==5) {legend('topright',legend=c('ERA5'),col=c('blue'),pch=15,cex=1.75,pt.cex=3)}
      if (i==6) {legend('topright',legend=c('PNWNAmet'),col=c('green'),pch=15,cex=1.75,pt.cex=3)}

      if (i==6) {
         axis(1,at=as.Date(c('1990-01-01','2000-01-01','2010-01-01','2018-01-01')),labels=c('1990','2000','2010','2018'),cex.axis=2.5,mgp=c(3,2,0))
      }
      box(which='plot',lwd=1.2)
    }

    ##Snowpack                   
    if (var.name=='pack') {                
      date.subset <- format(as.Date(dates),'%Y-%m-%d') %in% pillow.dates
      pillow.subset <- pillow.dates %in% format(as.Date(dates),'%Y-%m-%d')

      plot(as.Date(dates)[date.subset],apply(pack.sims[date.subset,],1,mean),type='l',lwd=3,col='white',cex.axis=1.5,yaxs='i',
               xlab='Date',ylab='Snowpack (cm)',main=site.names[i], cex.lab=1.75,cex.axis=1.75,cex.main=2)
      apply(pack.sims[date.subset,],2,function(x,y){lines(y,x,col='lightblue',lwd=2.5)},as.Date(dates[date.subset]))
      points(as.Date(pillow.dates),pillow.pack,cex=1,col='black',pch=16)
      lines(as.Date(dates[date.subset]),apply(pack.sims,1,mean)[date.subset],col='blue',lwd=3.5)
      abline(h=0)
    }
    }
}

mtext("Date",side=1,outer=TRUE,cex=2.0,line=4.6)
mtext("Snow Water Equivalent (mm)",side=2,outer=TRUE,cex=2.0,line=4.6)


    dev.off()