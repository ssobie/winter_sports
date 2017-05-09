##Script to compare snowdepth, SWE and snow cover between modelled values, 
##snow course sites and MODIS snow cover
library(ncdf4)

source('/storage/home/ssobie/code/repos/winter_sports/test.snow.model.r',chdir=T)

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


hyperbolic.snow <- function(pr,tasmax,tasmin,coeffs) {

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

##site <- 'callaghan'
model <- 'ERA'

b.coef.test <- seq(0.3,0.8,by=0.05)
c.coef.test <- seq(2.5,4.5,by=0.1)

b.min <- matrix(0,nrow=length(sites),ncol=12)
c.min <- matrix(0,nrow=length(sites),ncol=12)

for (s in seq_along(sites)) {
    site <- 'callaghan'
    print(s)
    ##Reanalysis 800m data
    clim.file <- paste0('/storage/data/projects/rci/data/winter_sports/BCCAQ2/snow_sites/',site,'_',model,'_800m_data.csv')
    clim.data <- read.csv(clim.file,header=T,as.is=T)
    pr.data <- clim.data$Pr
    tasmax.data <- clim.data$Tasmax
    tasmin.data <- clim.data$Tasmin
    tas.data <- clim.data$Tas
    dates <- clim.data$Dates

    ##Snow Course Data
    course.file <- paste('/storage/data/projects/rci/data/assessments/snow_model/',site,'_snow_course.csv',sep='')
    course.data <- read.csv(course.file,header=T,as.is=T)
    course.dates <- format(as.Date(course.data[,1]),'%Y-%m-%d')
    course.swe <- course.data[,3] ##mm
    course.pack <- course.data[,2] ##cm
    course.dense <-  course.data[,4]
    
    date.subset <- format(as.Date(dates),'%Y-%m-%d') %in% course.dates
    course.subset <- course.dates %in% format(as.Date(dates),'%Y-%m-%d') 
    
    coords <- get.coordinates(site)
    lat.bnds <- coords[2]
    elev <- coords[3]
    site.slope <- bc.slopes[which.min(abs(coords[1]-bc.lon)),which.min(abs(coords[2]-bc.lat))]    
    site.aspect <- bc.aspects[which.min(abs(coords[1]-bc.lon)),which.min(abs(coords[2]-bc.lat))]    


   ##---------------------------------------------
   for (z in 1:12) {

    if (1==1) {
    ##SWE Fitting
    dense.diff <- matrix(0,nrow=length(b.coef.test),ncol=length(c.coef.test)) 
    pack.diff <- matrix(0,nrow=length(b.coef.test),ncol=length(c.coef.test)) 
    swe.diff <- matrix(0,nrow=length(b.coef.test),ncol=length(c.coef.test)) 

    coeffs <- list(a=-49.49,b=0,c=0,d=1.0209)
    print(dim(swe.diff))
    for (i in seq_along(b.coef.test)) {
##        print(i)
        coeffs$b <- b.coef.test[i]
        for (j in seq_along(c.coef.test)) {
##          print(j)
          coeffs$c <- c.coef.test[j]        
          Snow <- hyperbolic.snow(pr.data,tasmax.data,tasmin.data,coeffs)
          results <- snow.melt(Date=dates, precip_mm=pr.data, Tmax_C=tasmax.data, Tmin_C=tasmin.data,Snow=Snow, 
                           lat_deg=lat.bnds, slope=site.slope, aspect=site.aspect, tempHt=1, windHt=2, groundAlbedo=0.25,
                           SurfEmissiv=0.95, windSp=1, forest=0, startingSnowDepth_m=0, startingSnowDensity_kg_m3=600)
          dense.diff[i,j] <- mean(course.dense[course.subset] - results$snowdense[date.subset],na.rm=T)
          pack.diff[i,j] <- mean(course.pack[course.subset] - results$snowdepth[date.subset]*1000,na.rm=T)
          swe.diff[i,j] <- mean(course.swe[course.subset] - results$swe[date.subset]*1000,na.rm=T)
      }
    }  
    }
 
    ix.min <- which(abs(swe.diff)==min(abs(swe.diff)),arr.ind=T)
    b.min[s,z] <- b.coef.test[ix.min[1]]
    c.min[s,z] <- c.coef.test[ix.min[2]]
    }
    }

    write.table(b.min,file='/storage/data/projects/rci/data/winter_sports/era_hypersnow_bcoef_fit.csv',sep=',',quote=F,row.names=sites)
    write.table(c.min,file='/storage/data/projects/rci/data/winter_sports/era_hypersnow_ccoef_fit.csv',sep=',',quote=F,row.names=sites)

    ##---------------------------------------------
    if (1==0) {
    ##Density Fitting
    DCoef.test <- seq(8,12,by=0.05)
    dense.diff <- rep(0,length(DCoef.test)) 
    pack.diff <- rep(0,length(DCoef.test)) 
    swe.diff <- rep(0,length(DCoef.test)) 

    for (i in seq_along(DCoef.test)) {
        results <- snow.melt(Date=dates, precip_mm=pr.data, Tmax_C=tasmax.data, Tmin_C=tasmin.data,DCoef.val=DCoef.test[i], 
                         lat_deg=lat.bnds, slope=site.slope, aspect=site.aspect, tempHt=1, windHt=2, groundAlbedo=0.25,
                         SurfEmissiv=0.95, windSp=1, forest=0, startingSnowDepth_m=0, startingSnowDensity_kg_m3=600)
        dense.diff[i] <- mean(course.dense[course.subset] - results$snowdense[date.subset],na.rm=T)
        pack.diff[i] <- mean(course.pack[course.subset] - results$snowdepth[date.subset]*100,na.rm=T)
        swe.diff[i] <- mean(course.swe[course.subset] - results$swe[date.subset]*1000,na.rm=T)

    }

    par(mfrow=c(2,1))    
    plot(DCoef.test,dense.diff)
    abline(v=DCoef.test[which.min(dense.diff)])
    abline(h=0)

    plot(DCoef.test,swe.diff)
    abline(v=DCoef.test[which.min(swe.diff)])
    abline(h=0)

    }


