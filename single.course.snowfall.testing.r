##Script to compare snowdepth, SWE and snow cover between modelled values, 
##snow course sites and MODIS snow cover
library(ncdf4)
library(plotrix)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)
source('/storage/home/ssobie/code/repos/winter_sports/site.coordinates.r',chdir=T)
source('/storage/home/ssobie/code/repos/winter_sports/snow.model.calibrated.r',chdir=T)

get_manual_parameters <- function(site) {
##
    ##coeffs <- list(a=-48.49,b=0.7628,c=3.0,d=1.0209) ##Best manual guess
    coeffs <- list(callaghan=list(a=-46.25,b=0.875,c=2.8,d=1.0209), ##Callaghan
                   grouse_mountain=list(a=-51.5,b=0.81,c=3.3,d=1.0209), ##Grouse Mountain
                   orchid_lake=list(a=-51.8,b=0.7,c=4.2,d=1.0209), ##Orchid Lake
                   palisade_lake=list(a=-51.85,b=0.35,c=4.45,d=1.0209), ##Palisade Lake                   
                   wahleach=list(a=-40.05,b=0.2,c=2.75,d=1.0209), ##Wahleach
                   brookmere=list(a=-40.0,b=0.05,c=2.7,d=1.0209), ##Brookmere
                   shovelnose_mountain=list(a=-40.15,b=0.075,c=2.75,d=1.0209), ##Shovelnose Mountain
                   lightning_lake=list(a=-40.15,b=0.05,c=2.7,d=1.0209), ##Lightning Lake
                   hamilton_hill=list(a=-40.45,b=0.1,c=2.85,d=1.0209), ##Hamilton Hill
                   klesilkwa=list(a=-40.1,b=0.175,c=2.67,d=1.0209), ##Klesilkwa
                   stave_lake=list(a=-50.6,b=0.825,c=3.4,d=1.0209), ##Stave Lake
                   dog_mountain=list(a=-50.5,b=0.725,c=3.03,d=1.0209), ##Dog Mountain
                   nahatlatch=list(a=-45.45,b=0.8,c=3.08,d=1.0209)) ##Nahatlatch

    rv <- coeffs[[site]]
    return(rv)      
}

get_cal_parameters <- function(site,model) {
   cal.file <- paste0('/storage/data/projects/rci/data/winter_sports/calibration_parameters/',site,'_',model,'_more_limited.csv')
   cal.data <- read.csv(cal.file,header=T,as.is=T)[[1]]
   coeffs <- list(a=cal.data[1],b=cal.data[2],c=cal.data[3],d=1.0209)
   return(coeffs)
}


get_parameters <- function(coords,type) {
   par.dir <- '/storage/data/projects/rci/data/winter_sports/'
   scale.nc <- nc_open(paste0(par.dir,'scale_hyper_snow_calibrated_',type,'.nc'))
   scales <- ncvar_get(scale.nc,'scale')
   lon <- ncvar_get(scale.nc,'lon')
   lat <- ncvar_get(scale.nc,'lat')
   nc_close(scale.nc)
   scale <- scales[which.min(abs(coords[1]-lon)),which.min(abs(coords[2]-lat))]

   slope.nc <- nc_open(paste0(par.dir,'slope_hyper_snow_calibrated_',type,'.nc'))
   slopes <- ncvar_get(slope.nc,'slope')
   slope <- slopes[which.min(abs(coords[1]-lon)),which.min(abs(coords[2]-lat))]
   nc_close(slope.nc)

   freq.nc <- nc_open(paste0(par.dir,'freq_hyper_snow_calibrated_',type,'.nc'))
   freqs <- ncvar_get(freq.nc,'freq')
   freq <- freqs[which.min(abs(coords[1]-lon)),which.min(abs(coords[2]-lat))]
   nc_close(freq.nc)

   coeffs <- list(a=-1*scale,b=slope,c=freq,d=1.0209)
   return(coeffs)
}

hyper_snow_test <- function(pr,tasmax,tasmin,coeffs) {

    ##coeffs <- list(a=-51.5,b=0.81,c=3.3,d=1.0209) ##Grouse Mountain    
    tas <- (tasmax+tasmin)/2                # degrees C
    frac <- coeffs$a*(tanh(coeffs$b*(tas-coeffs$c))-coeffs$d)
    sample <- runif(length(tas),min=0,max=100)
    test <- sample > frac
    high.temp <- tas > 15
    test[high.temp] <- TRUE
    low.temp <- tas < -10
    test[low.temp] <- FALSE

    snow.type <- rep(TRUE,length(tas))
    snow.type[test] <- FALSE

    NewSnowWatEq <- pr
    NewSnowWatEq[!snow.type] <- 0
    R_mm <- pr
    R_mm[snow.type] <- 0
    rv <- list(swe=NewSnowWatEq,
               rain=R_mm)
    return(rv)
}



get_snow_season_dates <- function(input.snow,dates) {

  years <- format(as.Date(dates),'%Y')
  jdays <- format(as.Date(dates),'%j')

  peaks.ix <- tapply(input.snow,as.factor(years),function(x){which.max(x[1:300])})
  peaks <- as.Date(paste0(levels(as.factor(years)),'-',sprintf('%03d',peaks.ix)),format='%Y-%j')
  
##Fix to return the date of peak SWE
  sim.snow <- input.snow
  sim.snow[sim.snow >= 25] <- 100
  sim.snow[sim.snow < 25] <- 0

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
  len <- length(ix)
  lengths <- ix[seq(2,len,2)] - ix[seq(1,len-1,2)]
  ###starts <- as.numeric(format(as.Date(dates[ix[seq(1,len-1,2)]]),'%j'))
  ###ends <- as.numeric(format(as.Date(dates[ix[seq(2,len,2)]]),'%j'))
  starts <- as.Date(dates[ix[seq(1,len-1,2)]])
  ends <- as.Date(dates[ix[seq(2,len,2)]])

  peaks.match <- as.Date('1900-01-01')
  for (k in seq_along(starts)) {
     peaks.match <- c(peaks.match,peaks[which.min(abs(peaks-starts[k]))])
  }

  rv <- list(lengths=lengths,
             starts =starts,
             peaks = peaks.match[-1],
             ends = ends,
             years = unique(years))
  return(rv)
}


read_snow_sim <- function(site,model,type) {

    ##PNWNAmet PRISM calibration
    model.dir <- '/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow_series/'
    swe.file <- paste0(model.dir,site,'_',model,'_',type,'_snow_model_data.csv')
    swe.data <- read.csv(swe.file,header=T,as.is=T)
    swe.values <- swe.data$SWE*1000
    swe.dates <- as.Date(swe.data$Dates)
    rv <- list(dates=swe.dates,swe=swe.values)
    return(rv)
}


##North Shore Sites

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

sites <- 'highland_valley'
inactive <- FALSE
course.dir <- 'snow_courses'

for (site in sites) {
site.name <- site

model <- 'PNWNAmet'

coords <- get_coordinates(site)
lat.bnds <- coords[2]
elev <- coords[3]


print(site)
##Reanalysis 800m data
clim.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/TPS/snow/snow_sites/',site,'_',model,'_800m_data.csv')
clim.data <- read.csv(clim.file,header=T,as.is=T)
clim.pr <- clim.data$Pr
clim.tasmax <- clim.data$Tasmax
clim.tasmin <- clim.data$Tasmin
clim.tas <- clim.data$Tas
clim.dates <- clim.data$Dates

Tmax.check <- clim.tasmax <  clim.tasmin
clim.tasmax[Tmax.check] <- clim.tasmin[Tmax.check]+1

##Apply snow fall model
##Elevation combined coefficients
slen <- 10
snow.series <- matrix(0,nrow=length(clim.pr),ncol=slen)
rain.series <- matrix(0,nrow=length(clim.pr),ncol=slen)
coeffs <- get_parameters(coords,type='parameter_PNWNAmet_prism')           
##coeffs <- get_point_parameters(site)           
for (k in 1:slen) {
   snow <- hyper_snow_test(clim.pr,clim.tasmax,clim.tasmin,coeffs)
   snow.series[,k] <- snow$swe
   rain.series[,k] <- snow$rain
}
snow.series.mean <- apply(snow.series,1,mean)/1000
rain.series.mean <- apply(rain.series,1,mean)/1000
rain.series.mean[snow.series.mean !=0] <- 0 

##----------------------------------------------------------
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

site.slope <- bc.slopes[which.min(abs(coords[1]-bc.lon)),which.min(abs(coords[2]-bc.lat))]
site.aspect <- bc.aspects[which.min(abs(coords[1]-bc.lon)),which.min(abs(coords[2]-bc.lat))]

coeff.matrix <- matrix(0,nrow=3,ncol=4)

if (inactive) {
  coeff.matrix[1,] <- rep(NA,4)
} else {
coeffs <- get_cal_parameters(site,model)
coeff.matrix[1,] <- unlist(coeffs)

cal.results <- snow_melt(Date=clim.dates, precip_mm=clim.pr, Tmax_C=clim.tasmax, Tmin_C=clim.tasmin,
                              cal_scale=coeffs$a,cal_slope=coeffs$b,cal_freq=coeffs$c,
                              lat_deg=coords[2], slope=site.slope, aspect=site.aspect, tempHt=1, windHt=2, groundAlbedo=0.25,
                              SurfEmissiv=0.95, windSp=1, forest=0, startingSnowDepth_m=0, startingSnowDensity_kg_m3=600)
}
##Kriged Parameters
if (1==0) {
coeffs <- get_parameters(coords,type=paste0('parameter_',model,'_prism'))           ##get_point_parameters(site) ##
coeff.matrix[2,] <- unlist(coeffs)           
kriged.results <- snow_melt(Date=clim.dates, precip_mm=clim.pr, Tmax_C=clim.tasmax, Tmin_C=clim.tasmin,
                                 cal_scale=coeffs$a,cal_slope=coeffs$b,cal_freq=coeffs$c,
                                 lat_deg=coords[2], slope=site.slope, aspect=site.aspect, tempHt=1, windHt=2, groundAlbedo=0.25,
                                 SurfEmissiv=0.95, windSp=1, forest=0, startingSnowDepth_m=0, startingSnowDensity_kg_m3=600)

##TPS
coeffs <- get_parameters(coords,type=paste0('parameter_',model,'_prism_TPS'))
coeff.matrix[3,] <- unlist(coeffs)           
tps.results <- snow_melt(Date=clim.dates, precip_mm=clim.pr, Tmax_C=clim.tasmax, Tmin_C=clim.tasmin,
                                 cal_scale=coeffs$a,cal_slope=coeffs$b,cal_freq=coeffs$c,
                                 lat_deg=coords[2], slope=site.slope, aspect=site.aspect, tempHt=1, windHt=2, groundAlbedo=0.25,
                                 SurfEmissiv=0.95, windSp=1, forest=0, startingSnowDepth_m=0, startingSnowDensity_kg_m3=600)
}

##------------------------------------------------------------------------


##Snow Course Data
course.file <- paste0('/storage/data/projects/rci/data/winter_sports/obs/',course.dir,'/',site,'_snow_course.csv',sep='')
##course.file <- paste('/storage/data/projects/rci/data/assessments/snow_model)
course.data <- read.csv(course.file,header=T,as.is=T)
course.dates <- format(as.Date(course.data[,1]),'%Y-%m-%d')
course.swe <- course.data[,3] ##mm
course.range <- range(course.dates)

if (inactive) {
   cal.results <- tps.results
}
season.dates <- get_snow_season_dates(cal.results$swe*1000,clim.dates)

se.len <- length(season.dates$lengths)
print(season.dates$starts)
cumulative.sim.pr <- clim.data$Pr*0
cumulative.end.pr <- clim.data$Pr*0
cumulative.snow.mean <- clim.data$Pr*0
cumulative.snow.end <- clim.data$Pr*0

for (j in 1:se.len) {
   sst <- grep(season.dates$starts[j],clim.dates)
   sen <- grep(season.dates$peaks[j],clim.dates)
   een <- grep(season.dates$ends[j],clim.dates)

   cumulative.sim.pr[sst:sen] <- cumsum(clim.data$Pr[sst:sen])
   cumulative.end.pr[sst:een] <- cumsum(clim.data$Pr[sst:een])

   cumulative.snow.mean[sst:sen] <- cumsum((cal.results$snowseries)[sst:sen])
   cumulative.snow.end[sst:een] <- cumsum((cal.results$snowseries)[sst:een])
}

##Snow Course Comparison
yupp <- max(c(max(course.swe,na.rm=T),max(cumulative.sim.pr)))
ymax <- max(c(max(course.swe,na.rm=T),max(cumulative.sim.pr)))

dst <- '1981-10-01'
den <- '1990-06-30'

par(mfrow=c(2,1))
par(mar=c(4.1,5,1.1,1.1))
if (1==0) {
plot(as.Date(course.dates),course.swe,cex=1.5,col='black',pch=16,
     xlim=c(as.Date(dst),as.Date(den)),ylim=c(0,ymax),
     main='',xlab='Date',ylab='SWE (mm)', cex.lab=1.25,axes=F)
axis(1,at=as.Date(c('1998-01-01','2000-01-01','2002-01-01','2004-01-01')),label=c('1998','2000','2002','2004'),cex.axis=1.25)
axis(2,at=seq(0,yupp,round((yupp-100)/3,-2)),label=seq(0,yupp,round((yupp-100)/3,-2)),cex.axis=1.25)

lines(as.Date(clim.dates),cumulative.snow.mean,col='blue',lwd=2)
lines(as.Date(clim.dates),cumulative.snow.end,col='blue',lwd=2,lty=2)

points(as.Date(course.dates),course.swe,cex=1.75,col='black',pch=16)
text(as.Date('2005-01-01'),0.95*ymax,site.name,cex=2)
legend('topright',legend=c('Course Obs.','ERA','PNWNAmet'),col=c('black','blue','green'),pch=16,cex=1.15)

box(which='plot')
abline(h=0)

plot(as.Date(clim.dates),clim.tas,cex=1.5,col='white',pch=16,
     xlim=c(as.Date(dst),as.Date(den)),ylim=range(c(min(clim.tasmin),max(clim.tasmax))),
     main='',xlab='Date',ylab='Temperature (degC)', cex.lab=1.25,axes=T)
axis(1,at=as.Date(c('1998-01-01','2000-01-01','2002-01-01','2004-01-01')),label=c('1998','2000','2002','2004'),cex.axis=1.25)
##axis(2,at=seq(0,yupp,round((yupp-100)/3,-2)),label=seq(0,yupp,round((yupp-100)/3,-2)),cex.axis=1.25)
lines(as.Date(clim.dates),clim.tasmax,col='red',lwd=2)
lines(as.Date(clim.dates),clim.tas,col='orange',lwd=2)
lines(as.Date(clim.dates),clim.tasmin,col='goldenrod',lwd=2)
box(which='plot')
abline(h=0)
}


plot(as.Date(course.dates),course.swe,cex=1.5,col='black',pch=16,
     xlim=as.Date(course.range), ##xlim=c(as.Date('1970-01-01'),as.Date('2020-12-31')),
     ylim=c(0,ymax),
     main='',xlab='Date',ylab='SWE (mm)', cex.lab=1.25,axes=F)
axis(1,at=as.Date(c('1960-01-01','1970-01-01','1980-01-01','1990-01-01','2000-01-01','2010-01-01')),
       label=c('1960','1970','1980','1990','2000','2010'),cex.axis=1.25)
axis(2,at=seq(0,yupp,round((yupp-100)/3,-2)),label=seq(0,yupp,round((yupp-100)/3,-2)),cex.axis=1.25)

points(as.Date(course.dates),course.swe,cex=1.75,col='black',pch=16)
text(as.Date('2005-01-01'),0.95*ymax,site.name,cex=2)
legend('topright',legend=c('Course Obs.','Cal','Kriged','TPS'),col=c('black','blue','green','goldenrod'),pch=16,cex=1.15)
lines(as.Date(clim.dates),cal.results$swe*1000,col='blue',lwd=2)
##lines(as.Date(clim.dates),kriged.results$swe*1000,col='green',lwd=2)
##lines(as.Date(clim.dates),tps.results$swe*1000,col='goldenrod',lwd=2)

box(which='plot')
abline(h=0)


paper.coeffs <- list(a=-49.49,b=0.5628,c=1.5,d=1.0209)
tas <- seq(-20,20,0.05)

paper.frac <- paper.coeffs$a*(tanh(paper.coeffs$b*(tas-paper.coeffs$c))-paper.coeffs$d)
frac2 <- coeffs$a*(tanh(coeffs$b*(tas-coeffs$c))-coeffs$d)

par(mar=c(5,5,2,1))
plot(tas,paper.frac,type='l',lwd=2,xlab="Daily Average Temperature (\u00B0C)",ylab='Percent Snow (%)',
              xaxs='i',yaxs='i',cex.lab=2,cex.axis=2)
lines(tas,frac2,col='red',lwd=2)
box(which='plot')



clim.subset <- as.Date(clim.dates) %in% as.Date(course.dates)
course.subset <- as.Date(course.dates) %in% as.Date(clim.dates)

print(paste0(model,' Cal MAE point fitted'))
print(mean((cal.results$swe[clim.subset]*1000 - course.swe[course.subset]),na.rm=T))

##print(paste0(model,' Kriged MAE fitted'))
##print(mean((kriged.results$swe[clim.subset]*1000 - course.swe[course.subset]),na.rm=T))

##print(paste0(model,' TPS MAE fitted'))
##print(mean((tps.results$swe[clim.subset]*1000 - course.swe[course.subset]),na.rm=T))

print('Coeff Matrix')
##row.names(coeff.matrix) <- c('Cal','Krige','TPS')
##print(coeff.matrix)


}